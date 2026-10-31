"""Vector-valued polynomial patches on a closed spatial surface mesh.

C0 and fixed vertices form an affine linear space. C1 means equality of the
*full* transverse spatial vectors after projection perpendicular to the
shared edge tangent; it is nonlinear because that tangent is also unknown.
All derivatives below come from SymPy polynomials, never finite differences.
"""
from dataclasses import dataclass
from functools import cached_property

import numpy as np
import scipy.linalg
from scipy.optimize import least_squares
import sympy as sp

from filmstar.polynomial_basis import monomials
from filmstar.planar.mesh import edge_t
from filmstar.spatial.mesh import SurfaceMesh

s, t = sp.symbols('s t')


def _nodes(count):
    return (1-np.cos(np.linspace(0, np.pi, count)))/2


@dataclass(frozen=True)
class _NormalOccurrence:
    first: np.ndarray
    second: np.ndarray
    sign: float


class SurfaceSpace:
    def __init__(self, mesh: SurfaceMesh, degree: int, smoothness: int):
        if not isinstance(mesh, SurfaceMesh):
            raise TypeError('mesh must be a spatial SurfaceMesh')
        if not isinstance(degree, int) or degree < 1 or smoothness not in (0, 1):
            raise ValueError('Use degree >= 1 and smoothness 0 or 1')
        self.mesh, self.degree, self.smoothness = mesh, degree, smoothness
        self.faces = tuple(mesh.parameter_mesh.faces)
        self.face_index = {face.id: k for k, face in enumerate(self.faces)}
        self.polynomials = monomials(degree, (s, t))
        self.derivatives_s = tuple(sp.diff(value, s) for value in self.polynomials)
        self.derivatives_t = tuple(sp.diff(value, t) for value in self.polynomials)
        self._functions = tuple(tuple(sp.lambdify((s, t), value, 'numpy')
                                      for value in expressions)
                                for expressions in (self.polynomials,
                                                    self.derivatives_s,
                                                    self.derivatives_t))
        self.columns = len(self.faces)*len(self.polynomials)
        self._assemble_c0()

    def _values(self, ss, tt, derivative=0):
        ss, tt = np.broadcast_arrays(np.asarray(ss, dtype=float),
                                     np.asarray(tt, dtype=float))
        return np.stack([np.broadcast_to(function(ss, tt), ss.shape).ravel()
                         for function in self._functions[derivative]], axis=1)

    def _face_row(self, fid, ss, tt, derivative=0):
        values = self._values(np.asarray([ss]), np.asarray([tt]), derivative)[0]
        row = np.zeros(self.columns)
        offset = self.face_index[fid]*len(self.polynomials)
        row[offset:offset+len(values)] = values
        return row

    def _edge_local(self, edge, side, parameters):
        expressions = edge.source_local if side == 'source' else edge.target_local
        functions = [sp.lambdify(edge_t, value, 'numpy') for value in expressions]
        return tuple(np.broadcast_to(function(parameters), parameters.shape).astype(float)
                     for function in functions)

    def _assemble_c0(self):
        rows = []
        targets = []
        nodes = _nodes(self.degree+1)
        for edge in self.mesh.edges:
            source = self._edge_local(edge, 'source', nodes)
            target = self._edge_local(edge, 'target', nodes)
            for k in range(len(nodes)):
                rows.append(self._face_row(edge.source_face, source[0][k], source[1][k])
                            - self._face_row(edge.target_face, target[0][k], target[1][k]))
                targets.append(np.zeros(3))
        anchors = []
        origins = {face.id: self.mesh.parameter_mesh.vertices[face.local_origin]
                   for face in self.faces}
        for vertex_id in self.mesh.vertex_ids:
            original = next(vid for vid, representative in self.mesh.vertex_classes.items()
                            if representative == vertex_id)
            face = next(face for face in self.faces if original in face.corners)
            point = self.mesh.parameter_mesh.vertices[original]
            origin = origins[face.id]
            rows.append(self._face_row(face.id, float(point[0]-origin[0]),
                                       float(point[1]-origin[1])))
            position = np.asarray(self.mesh.spatial_vertices[original], dtype=float)
            anchors.append(position)
            targets.append(position)
        self.edge_row_count = len(rows)-len(anchors)
        if self.smoothness == 1 and self.mesh.radius is not None:
            for edge in self.mesh.edges:
                local = edge.source_local
                local_derivatives = [sp.diff(value, edge_t) for value in local]
                for endpoint, vid in enumerate(edge.source_vertices):
                    point = [float(value.subs(edge_t, endpoint)) for value in local]
                    tangent = [float(value.subs(edge_t, endpoint))
                               for value in local_derivatives]
                    rows.append(tangent[0]*self._face_row(edge.source_face,*point,derivative=1)
                                + tangent[1]*self._face_row(edge.source_face,*point,derivative=2))
                    first = np.asarray(self.mesh.spatial_vertices[edge.source_vertices[0]],
                                       dtype=float)
                    second = np.asarray(self.mesh.spatial_vertices[edge.source_vertices[1]],
                                        dtype=float)
                    radial = self.mesh.vertex_normals[self.mesh.vertex_classes[vid]]
                    chord = second-first
                    targets.append(chord-np.dot(chord, radial)*radial)
        self.c0_matrix = np.stack(rows)
        self.anchor_positions = np.stack(anchors)
        U, singular, Vh = scipy.linalg.svd(self.c0_matrix, full_matrices=True)
        tolerance = max(self.c0_matrix.shape)*np.finfo(float).eps*singular[0]*10
        rank = int(np.sum(singular > tolerance))
        rhs = np.stack(targets)
        self.c0_targets = rhs
        self.particular = ((Vh[:rank].T / singular[:rank]) @
                           (U[:, :rank].T @ rhs)).T
        self.nullspace = Vh[rank:].T
        if np.max(np.abs(self.c0_matrix@self.particular.T-rhs)) > 1e-8:
            raise ValueError('Fixed spatial vertices are incompatible with C0 edges')

    @property
    def c0_dimension(self):
        return 3*self.nullspace.shape[1]

    def coefficients(self, parameters):
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (3, self.nullspace.shape[1]):
            raise ValueError(f'Expected parameter array shape (3,{self.nullspace.shape[1]})')
        return self.particular + parameters@self.nullspace.T

    def _affine_seed(self):
        coefficients = np.zeros((3, self.columns))
        count = len(self.polynomials)
        for face in self.faces:
            origin = np.asarray(self.mesh.parameter_mesh.vertices[face.local_origin], dtype=float)
            flat = np.asarray([self.mesh.parameter_mesh.vertices[v] for v in face.corners],
                              dtype=float)-origin
            spatial = np.asarray([self.mesh.spatial_vertices[v] for v in face.corners],
                                 dtype=float)
            affine = np.linalg.lstsq(np.column_stack((np.ones(len(flat)), flat)),
                                     spatial, rcond=None)[0]
            start = self.face_index[face.id]*count
            coefficients[:, start:start+3] = affine.T
        return coefficients

    def c0_surface(self):
        seed = self._affine_seed()
        parameters = (seed-self.particular)@self.nullspace
        return SurfaceSolution(self, self.coefficients(parameters), parameters)

    @cached_property
    def _c1_data(self):
        count = 3*self.degree-2  # degree 3(d-1) residual needs 3d-2 nodes
        nodes = _nodes(count)
        tangent_rows, transverse_rows = [], []
        for edge in self.mesh.edges:
            source = self._edge_local(edge, 'source', nodes)
            target = self._edge_local(edge, 'target', nodes)
            tangent_functions = [sp.lambdify(edge_t, sp.diff(value, edge_t), 'numpy')
                                 for value in edge.source_local]
            ds, dt = (np.broadcast_to(fn(nodes), nodes.shape)
                      for fn in tangent_functions)
            for k in range(len(nodes)):
                tangent_rows.append(
                    ds[k]*self._face_row(edge.source_face, source[0][k], source[1][k], 1)
                    + dt[k]*self._face_row(edge.source_face, source[0][k], source[1][k], 2))
                row = np.zeros(self.columns)
                for side, fid, local, inward in (
                        ('source', edge.source_face, source, edge.source_inward),
                        ('target', edge.target_face, target, edge.target_inward)):
                    row += inward[0]*self._face_row(fid, local[0][k], local[1][k], 1)
                    row += inward[1]*self._face_row(fid, local[0][k], local[1][k], 2)
                transverse_rows.append(row)
        T = np.stack(tangent_rows)
        W = np.stack(transverse_rows)
        vertex_occurrences = {}
        parameter_vertices = self.mesh.parameter_mesh.vertices
        for face in self.faces:
            corners = face.corners
            origin = np.asarray(parameter_vertices[face.local_origin], dtype=float)
            for k, vid in enumerate(corners):
                q = np.asarray(parameter_vertices[vid], dtype=float)
                directions = [np.asarray(parameter_vertices[corners[(k+step)%len(corners)]],
                                         dtype=float)-q for step in (1, -1)]
                directions = [direction/np.linalg.norm(direction)
                              for direction in directions]
                determinant = np.linalg.det(np.stack(directions))
                sign = self.mesh.face_orientation[face.id]*np.sign(determinant)
                local = q-origin
                first = (directions[0][0]*self._face_row(face.id,*local,derivative=1)
                         + directions[0][1]*self._face_row(face.id,*local,derivative=2))
                second = (directions[1][0]*self._face_row(face.id,*local,derivative=1)
                          + directions[1][1]*self._face_row(face.id,*local,derivative=2))
                vertex_occurrences.setdefault(self.mesh.vertex_classes[vid], []).append(
                    _NormalOccurrence(first, second, sign))
        normal_pairs = [(occurrences[0], occurrence)
                        for occurrences in vertex_occurrences.values()
                        for occurrence in occurrences[1:]]
        # Face samples choose an interior grid; the surface-radius objective
        # helps select a visually useful member of the nonlinear solution set.
        radius_rows, derivative_rows_s, derivative_rows_t, sample_signs = [], [], [], []
        for face in self.faces:
            points = np.asarray([parameter_vertices[v] for v in face.corners], dtype=float)
            origin = np.asarray(parameter_vertices[face.local_origin], dtype=float)
            centroid = np.mean(points, axis=0)
            for vertex in points:
                point = 0.75*centroid+0.25*vertex-origin
                radius_rows.append(self._face_row(face.id, point[0], point[1]))
                derivative_rows_s.append(self._face_row(face.id, point[0], point[1], 1))
                derivative_rows_t.append(self._face_row(face.id, point[0], point[1], 2))
                sample_signs.append(self.mesh.face_orientation[face.id])
        return (T, W, normal_pairs, np.stack(radius_rows),
                np.stack(derivative_rows_s), np.stack(derivative_rows_t),
                np.asarray(sample_signs))

    def c1_errors(self, coefficients):
        Tmat, Wmat, pairs, _, _, _, _ = self._c1_data
        tangents = (Tmat@coefficients.T)
        transverse = (Wmat@coefficients.T)
        length_squared = np.sum(tangents*tangents, axis=1)
        if np.min(length_squared) <= 1e-16:
            return float('inf'), float('inf'), 0.0
        projected = transverse-(np.sum(transverse*tangents, axis=1)/length_squared)[:, None]*tangents
        normal_error = 0.0
        min_normal = float('inf')
        for a, b in pairs:
            first = a.sign*np.cross(a.first@coefficients.T, a.second@coefficients.T)
            second = b.sign*np.cross(b.first@coefficients.T, b.second@coefficients.T)
            na, nb = np.linalg.norm(first), np.linalg.norm(second)
            min_normal = min(min_normal, na, nb)
            if min(na, nb) <= 1e-12:
                return float('inf'), float('inf'), 0.0
            normal_error = max(normal_error, np.linalg.norm(first/na-second/nb))
        return float(np.max(np.linalg.norm(projected, axis=1))), normal_error, min_normal

    def solve_c1(self, *, max_nfev=120, verbose=0, shape_weight=0.05,
                 perturb_scale=0.1):
        if self.smoothness != 1:
            raise ValueError('C1 solving needs smoothness=1')
        Tmat, Wmat, pairs, radius_matrix, Umat, Vmat, sample_signs = self._c1_data
        nullspace = self.nullspace
        free = nullspace.shape[1]
        Tz, Wz, Xz, Uz, Vz = (matrix@nullspace
                              for matrix in (Tmat, Wmat, radius_matrix, Umat, Vmat))
        seed = self._affine_seed()
        initial = (seed-self.particular)@nullspace
        # The affine seed can project to z=0 exactly. TRF then starts with an
        # almost zero trust radius and falsely reports convergence.
        choices = np.asarray([-1, -2/3, -1/2, -1/3, 1/3, 1/2, 2/3, 1])
        initial += perturb_scale*np.resize(choices, initial.shape)
        center = self.mesh.center
        radius = self.mesh.radius or float(np.median(
            [np.linalg.norm(np.asarray(point, dtype=float)-center)
             for point in self.mesh.spatial_vertices.values()]))
        reference_points = radius_matrix@seed.T
        reference_u, reference_v = Umat@seed.T, Vmat@seed.T
        reference_normal = sample_signs[:, None]*np.cross(reference_u, reference_v)
        reference_area = np.linalg.norm(reference_normal, axis=1)
        reference_unit = reference_normal/reference_area[:, None]
        minimum_area = 0.2*reference_area
        sphere_weight = shape_weight

        def evaluate(flat, with_jacobian=False):
            parameters = flat.reshape(3, free)
            coefficients = self.coefficients(parameters)
            tangents = Tmat@coefficients.T
            transverse = Wmat@coefficients.T
            length2 = np.maximum(np.sum(tangents*tangents, axis=1), 1e-12)
            dot = np.sum(tangents*transverse, axis=1)
            alpha = dot/length2
            edge_residual = transverse-alpha[:, None]*tangents
            normals = []
            normal_data = []
            for a, b in pairs:
                pair_normals = []
                pair_data = []
                for occurrence in (a, b):
                    p = occurrence.first@coefficients.T
                    q = occurrence.second@coefficients.T
                    raw = occurrence.sign*np.cross(p, q)
                    length = max(np.linalg.norm(raw), 1e-12)
                    pair_normals.append(raw/length)
                    pair_data.append((p, q, raw/length, length, occurrence))
                normals.append(pair_normals[1]-pair_normals[0])
                normal_data.append(pair_data)
            raw_points = radius_matrix@coefficients.T
            points = raw_points-center
            distances = np.maximum(np.linalg.norm(points, axis=1), 1e-12)
            if self.mesh.radius is None:
                sphere = sphere_weight*(raw_points-reference_points).ravel()/radius
            else:
                sphere = sphere_weight*(distances-radius)/radius
            partial_u, partial_v = Umat@coefficients.T, Vmat@coefficients.T
            normal = sample_signs[:, None]*np.cross(partial_u, partial_v)
            signed_area = np.sum(normal*reference_unit, axis=1)
            barrier = np.minimum(signed_area-minimum_area, 0)/minimum_area
            residual = np.concatenate((edge_residual.ravel(),
                                       np.asarray(normals).ravel(), sphere, barrier))
            if not with_jacobian:
                return residual
            jacobian = np.zeros((len(residual), 3*free))
            edge_count = len(tangents)
            for axis in range(3):
                dT, dW = Tz, Wz
                d_dot = dT*transverse[:, axis, None] + tangents[:, axis, None]*dW
                d_length2 = 2*tangents[:, axis, None]*dT
                d_alpha = (d_dot-alpha[:, None]*d_length2)/length2[:, None]
                block = -tangents[:, :, None]*d_alpha[:, None, :]
                block[:, axis, :] += dW-alpha[:, None]*dT
                jacobian[:3*edge_count, axis*free:(axis+1)*free] = block.reshape(-1, free)
                for k, (a_data, b_data) in enumerate(normal_data):
                    variation = []
                    for p, q, unit, length, occurrence in (a_data, b_data):
                        basis = np.eye(3)[axis]
                        dp = occurrence.first@nullspace
                        dq = occurrence.second@nullspace
                        d_raw = occurrence.sign*(
                            np.cross(dp[:, None]*basis, q)
                            + np.cross(p, dq[:, None]*basis))
                        d_unit = (d_raw-unit[None, :]*np.sum(d_raw*unit, axis=1)[:, None])/length
                        variation.append(d_unit)
                    start = 3*edge_count+3*k
                    jacobian[start:start+3, axis*free:(axis+1)*free] = (variation[1]-variation[0]).T
                start = 3*edge_count+3*len(normal_data)
                if self.mesh.radius is None:
                    block = np.zeros((len(reference_points), 3, free))
                    block[:, axis, :] = sphere_weight/radius*Xz
                    jacobian[start:start+len(sphere), axis*free:(axis+1)*free] = block.reshape(-1, free)
                else:
                    jacobian[start:start+len(sphere), axis*free:(axis+1)*free] = (
                        sphere_weight/radius*(points[:, axis]/distances)[:, None]*Xz)
                basis = np.eye(3)[axis]
                d_normal = sample_signs[:, None, None]*(
                    np.cross(Uz[:, :, None]*basis, partial_v[:, None, :])
                    + np.cross(partial_u[:, None, :], Vz[:, :, None]*basis))
                d_area = np.sum(d_normal*reference_unit[:, None, :], axis=-1)
                jacobian[start+len(sphere):, axis*free:(axis+1)*free] = (
                    (signed_area < minimum_area)[:, None]*d_area/minimum_area[:, None])
            return jacobian

        first = least_squares(lambda z: evaluate(z), initial.ravel(),
                               jac=lambda z: evaluate(z, True), method='lm',
                               max_nfev=max(30, max_nfev//3),
                               ftol=1e-11, xtol=1e-11, gtol=1e-11,
                               verbose=verbose)
        sphere_weight = 0.0
        result = least_squares(lambda z: evaluate(z), first.x,
                               jac=lambda z: evaluate(z, True), method='lm',
                               max_nfev=max(50, max_nfev-first.nfev),
                               ftol=1e-12, xtol=1e-12, gtol=1e-12,
                               verbose=verbose)
        parameters = result.x.reshape(3, free)
        for _ in range(30):
            current = self.coefficients(parameters)
            tangents = Tmat@current.T
            transverse = Wmat@current.T
            squared = np.sum(tangents*tangents, axis=1)
            if np.min(squared) < 1e-10:
                break
            shear = np.sum(tangents*transverse, axis=1)/squared
            linear = Wmat-shear[:, None]*Tmat
            restricted = linear@nullspace
            delta = np.stack([scipy.linalg.lstsq(restricted, -linear@current[axis],
                                                 cond=1e-10)[0]
                              for axis in range(3)])
            before = self.c1_errors(current)[0]
            for factor in (1.0, 0.5, 0.25, 0.125):
                candidate = parameters+factor*delta
                if self.c1_errors(self.coefficients(candidate))[0] < before:
                    parameters = candidate
                    break
            else:
                break
            if self.c1_errors(self.coefficients(parameters))[0] < 1e-7:
                break
        solution = SurfaceSolution(self, self.coefficients(parameters), parameters)
        edge_error, normal_error, min_normal = self.c1_errors(solution.coefficients)
        if edge_error > 1e-6 or normal_error > 1e-5 or min_normal < 1e-8:
            raise ValueError(f'C1 solve failed: edge={edge_error:.3g}, '
                             f'vertex={normal_error:.3g}, normal={min_normal:.3g}; '
                             f'{result.message}')
        if solution.independent_edge_error() > 1e-6:
            raise ValueError('C1 solve failed at independent edge samples')
        metrics = solution.geometry_metrics()
        if (metrics['minimum_normal_ratio'] < 1e-3
                or (metrics['minimum_outward_cosine'] is not None
                    and metrics['minimum_outward_cosine'] < 0.05)
                or (metrics['maximum_radius_deviation'] is not None
                    and metrics['maximum_radius_deviation'] > 0.3)):
            raise ValueError(f'C1 surface folds or collapses: {metrics}')
        return solution


@dataclass(frozen=True)
class SurfaceSolution:
    space: SurfaceSpace
    coefficients: np.ndarray
    parameters: np.ndarray

    def evaluate_local(self, face_id, ss, tt):
        start = self.space.face_index[face_id]*len(self.space.polynomials)
        face_coefficients = self.coefficients[:, start:start+len(self.space.polynomials)]
        return (self.space._values(ss, tt)@face_coefficients.T).reshape(
            np.broadcast_shapes(np.shape(ss), np.shape(tt))+(3,))

    def evaluate_local_derivatives(self, face_id, ss, tt):
        shape = np.broadcast_shapes(np.shape(ss), np.shape(tt))+(3,)
        start = self.space.face_index[face_id]*len(self.space.polynomials)
        face_coefficients = self.coefficients[:, start:start+len(self.space.polynomials)]
        return tuple((self.space._values(ss, tt, derivative)@face_coefficients.T).reshape(shape)
                     for derivative in (0, 1, 2))

    def geometry_metrics(self, samples=11):
        from filmstar.plotting.data import PolygonGrid

        lengths, facing, deviations = [], [], []
        mesh = self.space.mesh
        for face in self.space.faces:
            grid = PolygonGrid(face.id, tuple(mesh.parameter_mesh.vertices[v]
                                              for v in face.corners))
            xx, yy = grid.sample(np.linspace(0, 1, samples))
            origin = mesh.parameter_mesh.vertices[face.local_origin]
            xyz, first, second = self.evaluate_local_derivatives(
                face.id, xx-float(origin[0]), yy-float(origin[1]))
            normal = mesh.face_orientation[face.id]*np.cross(first, second)
            length = np.linalg.norm(normal, axis=-1)
            lengths.extend(length.ravel())
            if mesh.radius is not None:
                radial = xyz-mesh.center
                radial_length = np.linalg.norm(radial, axis=-1)
                facing.extend((np.sum(normal*radial, axis=-1)/
                               np.maximum(length*radial_length, 1e-12)).ravel())
                deviations.extend((np.abs(radial_length-mesh.radius)/mesh.radius).ravel())
        lengths = np.asarray(lengths)
        return {'minimum_normal_ratio': float(np.min(lengths)/np.median(lengths)),
                'minimum_outward_cosine': float(np.min(facing)) if facing else None,
                'maximum_radius_deviation': float(np.max(deviations)) if deviations else None}

    def independent_edge_error(self, samples=31):
        levels = np.linspace(0, 1, samples)
        largest = 0.0
        for edge in self.space.mesh.edges:
            source = self.space._edge_local(edge, 'source', levels)
            target = self.space._edge_local(edge, 'target', levels)
            _, su, sv = self.evaluate_local_derivatives(edge.source_face, *source)
            _, tu, tv = self.evaluate_local_derivatives(edge.target_face, *target)
            derivatives = [np.broadcast_to(
                sp.lambdify(edge_t, sp.diff(value, edge_t), 'numpy')(levels), levels.shape)
                for value in edge.source_local]
            tangent = (derivatives[0][:, None]*su
                       + derivatives[1][:, None]*sv)
            transverse = (edge.source_inward[0]*su+edge.source_inward[1]*sv
                          + edge.target_inward[0]*tu+edge.target_inward[1]*tv)
            squared = np.sum(tangent*tangent, axis=-1)
            if np.min(squared) < 1e-12:
                return float('inf')
            projected = transverse-(np.sum(tangent*transverse, axis=-1)/squared)[:, None]*tangent
            largest = max(largest, float(np.max(np.linalg.norm(projected, axis=-1))))
        return largest

    def evaluate_face(self, face_id, xx, yy):
        face = self.space.faces[self.space.face_index[face_id]]
        origin = self.space.mesh.parameter_mesh.vertices[face.local_origin]
        return self.evaluate_local(face_id, np.asarray(xx)-float(origin[0]),
                                   np.asarray(yy)-float(origin[1]))
