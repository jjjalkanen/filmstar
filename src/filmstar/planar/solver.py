"""Exact fixed-degree spline spaces on topology from :mod:`spline_mesh`.

Across an interior edge b=0, adjacent pieces differ by a multiple of
b**(smoothness+1). The symbolic kernel is the degree part of Lemma 2.1 in
the supplied DiPasquale--Sottile paper. Plot geometry is not used here.
"""
from dataclasses import dataclass
from functools import cached_property
import signal
import threading

import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix

from filmstar.polynomial_basis import monomials
from filmstar.planar.mesh import Mesh, x, y

s, t, u, v = sp.symbols('s t u v')


@dataclass(frozen=True)
class MatchingConstraint:
    id: str
    source_face: str
    target_face: str
    boundary: sp.Expr
    matrix: sp.ImmutableMatrix
    offset: tuple[sp.Expr, sp.Expr]
    anchor: tuple[sp.Expr, sp.Expr]


def _vertex_normal_rows(mesh, basis, columns):
    """Compare normals from two face tangents at each identified corner.

    For a regular height graph the normalized cross product of two lifted
    nonparallel edge tangents is (-p_x,-p_y,1). Rotating its horizontal part
    into a common face chart makes the equations linear in spline coefficients.
    """
    from collections import defaultdict, deque

    adjacency = defaultdict(list)
    identity = sp.ImmutableMatrix.eye(2)

    def connect(a, b, rotation):
        adjacency[a].append((b, rotation))
        adjacency[b].append((a, rotation.inv()))

    def matching_vertices(source, target, rotation, offset):
        for vertex in source.vertices:
            mapped = rotation*sp.Matrix(mesh.vertices[vertex])+sp.Matrix(offset)
            for other in target.vertices:
                if all(sp.simplify(mapped[k]-mesh.vertices[other][k]) == 0
                       for k in range(2)):
                    yield vertex, other
                    break
            else:
                raise ValueError(f'Could not pair vertices of {source.id}, {target.id}')

    for edge in mesh.interior_edges:
        for vertex in edge.vertices:
            connect((edge.faces[0], vertex), (edge.faces[1], vertex), identity)
    edges = {edge.id: edge for edge in mesh.edges}
    for gluing in mesh.gluings:
        source, target = edges[gluing.source], edges[gluing.target]
        for a, b in matching_vertices(source, target, gluing.matrix, gluing.offset):
            connect((source.faces[0], a), (target.faces[0], b), gluing.matrix)
    for seam in mesh.periodic_identifications:
        source, target = edges[seam.source], edges[seam.target]
        for a, b in matching_vertices(source, target, identity, seam.translation):
            connect((source.faces[0], a), (target.faces[0], b), identity)

    face_index = {face.id: k for k, face in enumerate(mesh.faces)}
    face_origin = {face.id: mesh.vertices[face.local_origin] for face in mesh.faces}
    faces = {face.id: face for face in mesh.faces}
    normal_parts = {}

    def face_normal_parts(fid, vid):
        """Cross two distinct boundary tangents, with vertical component one."""
        key = (fid, vid)
        if key in normal_parts:
            return normal_parts[key]
        corners = faces[fid].corners
        if not corners:
            raise ValueError(f'Face {fid} needs ordered corners for vertex normals')
        position = corners.index(vid)
        point = sp.Matrix(mesh.vertices[vid])
        directions = [sp.Matrix(mesh.vertices[corners[(position+step)%len(corners)]])-point
                      for step in (1, -1)]
        determinant = sp.expand(directions[0][0]*directions[1][1]
                                - directions[0][1]*directions[1][0])
        if determinant == 0:
            raise ValueError(f'Collinear incident edges at {fid}, {vid}')
        ox, oy = face_origin[fid]
        parts = []
        for monomial in basis:
            derivative = sp.Matrix([sp.diff(monomial, s).subs({s: point[0]-ox,
                                                                t: point[1]-oy}),
                                    sp.diff(monomial, t).subs({s: point[0]-ox,
                                                                t: point[1]-oy})])
            lifted = [sp.Matrix([direction[0], direction[1],
                                 derivative.dot(direction)])
                      for direction in directions]
            normal = lifted[0].cross(lifted[1])/determinant
            parts.append(sp.Matrix(normal[:2]))
        normal_parts[key] = parts
        return parts
    rows = []
    visited = set()
    for root in adjacency:
        if root in visited:
            continue
        paths = {root: identity}
        queue = deque([root])
        visited.add(root)
        while queue:
            node = queue.popleft()
            for neighbor, rotation in adjacency[node]:
                if neighbor not in paths:
                    paths[neighbor] = rotation*paths[node]
                    visited.add(neighbor)
                    queue.append(neighbor)
        if len(paths) < 2:
            continue
        for node in list(paths)[1:]:
            for component in range(2):
                row = [sp.S.Zero]*columns
                for sign, occurrence, rotation in ((1, node, paths[node]),
                                                   (-1, root, identity)):
                    fid, vid = occurrence
                    for k, normal in enumerate(face_normal_parts(fid, vid)):
                        row[face_index[fid]*len(basis)+k] += sign*sp.expand(
                            (rotation.T*normal)[component])
                rows.append(row)
    return sp.Matrix(rows) if rows else sp.zeros(0, columns)


class SplineSpace:
    def __init__(self, mesh: Mesh, degree: int, smoothness: int,
                 vertex_normals: bool | None = None):
        if not isinstance(mesh, Mesh):
            raise TypeError('mesh must be a loaded Mesh')
        if (not isinstance(degree, int) or not isinstance(smoothness, int)
                or degree < 0 or smoothness < 0):
            raise ValueError('degree and smoothness must be nonnegative integers')
        if vertex_normals is None:
            vertex_normals = bool(mesh.gluings and smoothness >= 1)
        self.mesh, self.degree, self.smoothness = mesh, degree, smoothness
        self.vertex_normals = vertex_normals
        self.monomials = monomials(degree, (s, t))
        by_id = {edge.id: edge for edge in mesh.edges}
        identity = sp.ImmutableMatrix.eye(2)
        self.constraints = tuple(
            MatchingConstraint(edge.id, edge.faces[1], edge.faces[0],
                               edge.polynomial, identity, (sp.S.Zero, sp.S.Zero),
                               mesh.vertices[edge.vertices[0]])
            for edge in mesh.interior_edges) + tuple(
            MatchingConstraint(item.id, by_id[item.source].faces[0],
                               by_id[item.target].faces[0],
                               by_id[item.source].polynomial, identity, item.translation,
                               mesh.vertices[by_id[item.source].vertices[0]])
            for item in mesh.periodic_identifications) + tuple(
            MatchingConstraint(item.id, by_id[item.source].faces[0],
                               by_id[item.target].faces[0],
                               by_id[item.source].polynomial, item.matrix, item.offset,
                               mesh.vertices[by_id[item.source].vertices[0]])
            for item in mesh.gluings)
        self.cofactor_monomials = tuple(
            monomials(degree - (smoothness+1)*sp.Poly(item.boundary, x, y).total_degree(),
                      (u, v))
            for item in self.constraints)
        n, nf, nc = len(self.monomials), len(mesh.faces), len(self.constraints)
        self.matrix = sp.zeros(nc*n, nf*n + sum(map(len, self.cofactor_monomials)))
        indices = {sp.Poly(m, s, t).monoms()[0]: k for k, m in enumerate(self.monomials)}
        face_index = {face.id: i for i, face in enumerate(mesh.faces)}
        face_origin = {face.id: mesh.vertices[face.local_origin] for face in mesh.faces}
        offset = nf*n
        for e, (item, qm) in enumerate(zip(self.constraints, self.cofactor_monomials)):
            source = face_index[item.source_face]
            target = face_index[item.target_face]
            mapped_anchor = item.matrix*sp.Matrix(item.anchor) + sp.Matrix(item.offset)
            ax, ay = item.anchor
            source_origin = face_origin[item.source_face]
            target_origin = face_origin[item.target_face]
            for k, m in enumerate(self.monomials):
                source_local = sp.Poly(m.subs(
                    {s: ax+u-source_origin[0], t: ay+v-source_origin[1]},
                    simultaneous=True), u, v)
                for powers, coefficient in source_local.terms():
                    self.matrix[e*n+indices[powers], source*n+k] -= coefficient
                target_local = sp.Poly(m.subs(
                    {s: mapped_anchor[0]+item.matrix[0, 0]*u+item.matrix[0, 1]*v-target_origin[0],
                     t: mapped_anchor[1]+item.matrix[1, 0]*u+item.matrix[1, 1]*v-target_origin[1]},
                    simultaneous=True), u, v)
                for powers, coefficient in target_local.terms():
                    self.matrix[e*n+indices[powers], target*n+k] += coefficient
            boundary_local = sp.expand(item.boundary.subs(
                {x: ax+u, y: ay+v}, simultaneous=True))
            for k, m in enumerate(qm):
                for powers, coefficient in sp.Poly(
                        boundary_local**(smoothness+1)*m, u, v).terms():
                    self.matrix[e*n+indices[powers], offset+k] = -coefficient
            offset += len(qm)
        self.edge_matrix = self.matrix
        if vertex_normals and smoothness >= 1:
            self.vertex_matrix = _vertex_normal_rows(mesh, self.monomials,
                                                     self.matrix.cols)
            self.matrix = self.matrix.col_join(self.vertex_matrix)
        else:
            self.vertex_matrix = sp.zeros(0, self.matrix.cols)

    @cached_property
    def kernel(self):
        if self.mesh.gluings and threading.current_thread() is threading.main_thread():
            previous_handler = signal.getsignal(signal.SIGALRM)
            def expired(_signum, _frame):
                raise TimeoutError('Exact nullspace exceeded 55 seconds')
            signal.signal(signal.SIGALRM, expired)
            previous_timer = signal.setitimer(signal.ITIMER_REAL, 55)
            try:
                return self._compute_kernel()
            finally:
                signal.setitimer(signal.ITIMER_REAL, *previous_timer)
                signal.signal(signal.SIGALRM, previous_handler)
        return self._compute_kernel()

    def _compute_kernel(self):
        radicals = sorted({power for point in self.mesh.vertices.values()
                           for coordinate in point for power in coordinate.atoms(sp.Pow)
                           if power.exp.is_Rational and power.exp.q > 1}, key=str)
        if radicals:
            field = sp.QQ.algebraic_field(*radicals)
            # Expression-domain elimination expands nested radicals at every pivot.
            # Work in their exact number field instead.
            converted = {}
            def to_domain(matrix):
                rows = {}
                for (i, j), value in matrix.todok().items():
                    if value not in converted:
                        converted[value] = field.from_sympy(value)
                    rows.setdefault(i, {})[j] = converted[value]
                return DomainMatrix(rows, matrix.shape, field)
            edge_domain = to_domain(self.edge_matrix)
            vectors = edge_domain.nullspace()
            domain_matrix = edge_domain
            if self.vertex_matrix.rows and vectors.shape[0]:
                # A degree-one C1 edge kernel can contain only constants.
                # Vertex-normal rows differentiate the face pieces, so in
                # that case their restriction to the edge kernel is exactly
                # zero. Test sparse support before converting algebraic
                # coefficients into the number field; that conversion can
                # dominate the whole solve for large polygonal cutouts.
                active_columns = {
                    column for row in vectors.to_list()
                    for column, value in enumerate(row) if value != field.zero
                }
                relevant = any(column in active_columns
                               for _, column in self.vertex_matrix.todok())
                if relevant:
                    vertex_domain = to_domain(self.vertex_matrix)
                    projected = vertex_domain*vectors.transpose()
                    combinations = projected.nullspace()
                    vectors = combinations*vectors
                    domain_matrix = edge_domain.vstack(vertex_domain)
            if vectors.shape[0] == 0:
                raise ValueError('Matching constraints admit no nonzero spline')
            face_columns = len(self.mesh.faces)*len(self.monomials)
            rows = vectors.to_list()
            normalized = []
            for row in rows:
                pivot = next((value for value in row[:face_columns]
                              if value != field.zero), None)
                if pivot is None:
                    raise ValueError('Null vector has no face polynomial')
                normalized.append([value/pivot for value in row])
            self._exact_field = field
            self._domain_matrix = domain_matrix
            return sp.ImmutableMatrix(DomainMatrix(
                normalized, vectors.shape, field).to_Matrix().T)
        vectors = self.matrix.nullspace()
        if not vectors:
            raise ValueError('Matching constraints admit no nonzero spline')
        return sp.ImmutableMatrix.hstack(*vectors)

    @property
    def dimension(self):
        return self.kernel.cols

    @property
    def rank(self):
        return self.matrix.cols - self.dimension

    def unpack(self, vector):
        n = len(self.monomials)
        pieces = {face.id: sp.expand(sum(vector[i*n+k]*m for k, m in enumerate(self.monomials)))
                  for i, face in enumerate(self.mesh.faces)}
        offset, cofactors = len(self.mesh.faces)*n, {}
        for item, qm in zip(self.constraints, self.cofactor_monomials):
            cofactors[item.id] = sp.expand(sum(vector[offset+k]*m for k, m in enumerate(qm)))
            offset += len(qm)
        return pieces, cofactors

    @cached_property
    def basis(self):
        return tuple(self.unpack(self.kernel[:, k])[0] for k in range(self.dimension))

    def spline(self, parameters):
        parameters = tuple(sp.Rational(value) for value in parameters)
        if len(parameters) != self.dimension:
            raise ValueError(f'Expected {self.dimension} basis coefficients')
        return Spline(self, parameters)


class Spline:
    def __init__(self, space: SplineSpace, parameters):
        self.space, self.parameters = space, parameters
        self._vector = space.kernel*sp.Matrix(parameters)
        self.pieces, self.cofactors = space.unpack(self._vector)
        self._functions = {fid: sp.lambdify((s, t), polynomial, 'numpy')
                           for fid, polynomial in self.pieces.items()}

    def evaluate_face(self, face_id, xx, yy):
        """Evaluate using physical coordinates; subtract the face corner first."""
        if face_id not in self.pieces:
            raise ValueError(f'Unknown face: {face_id}')
        xx, yy = np.broadcast_arrays(np.asarray(xx, dtype=float), np.asarray(yy, dtype=float))
        face = next(face for face in self.space.mesh.faces if face.id == face_id)
        ox, oy = self.space.mesh.vertices[face.local_origin]
        return self.evaluate_local(face_id, xx-float(ox), yy-float(oy))

    def evaluate_local(self, face_id, ss, tt):
        """Evaluate the face polynomial directly in its own affine s,t chart."""
        if face_id not in self.pieces:
            raise ValueError(f'Unknown face: {face_id}')
        ss, tt = np.broadcast_arrays(np.asarray(ss, dtype=float), np.asarray(tt, dtype=float))
        return np.broadcast_to(self._functions[face_id](ss, tt), ss.shape)

    def physical_polynomial(self, face_id):
        """Exact global expression for algebraic checks; numerical evaluation uses s,t."""
        face = next(face for face in self.space.mesh.faces if face.id == face_id)
        ox, oy = self.space.mesh.vertices[face.local_origin]
        return sp.expand(self.pieces[face_id].subs(
            {s: x-ox, t: y-oy}, simultaneous=True))

    def verify(self):
        if hasattr(self.space, '_domain_matrix'):
            field = self.space._exact_field
            vector = DomainMatrix.from_Matrix(self._vector).convert_to(field)
            return (self.space._domain_matrix*vector).is_zero_matrix
        def identity(item):
            ax, ay = item.anchor
            mapped_anchor = item.matrix*sp.Matrix(item.anchor) + sp.Matrix(item.offset)
            origins = {face.id: self.space.mesh.vertices[face.local_origin]
                       for face in self.space.mesh.faces}
            so, to = origins[item.source_face], origins[item.target_face]
            source = self.pieces[item.source_face].subs(
                {s: ax+u-so[0], t: ay+v-so[1]}, simultaneous=True)
            target = self.pieces[item.target_face].subs(
                {s: mapped_anchor[0]+item.matrix[0, 0]*u+item.matrix[0, 1]*v-to[0],
                 t: mapped_anchor[1]+item.matrix[1, 0]*u+item.matrix[1, 1]*v-to[1]},
                simultaneous=True)
            boundary = item.boundary.subs({x: ax+u, y: ay+v}, simultaneous=True)
            return sp.expand(target-source
                             - boundary**(self.space.smoothness+1)*self.cofactors[item.id])
        return all(identity(item) == 0 for item in self.space.constraints)
