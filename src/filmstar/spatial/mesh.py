"""Spatial surface data over a separate, unchanged planar parameter mesh."""
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
import sympy as sp

from filmstar.planar.mesh import Mesh, edge_t, expression, load_mesh


@dataclass(frozen=True)
class SurfaceEdge:
    id: str
    source_vertices: tuple[str, str]
    source_face: str
    target_face: str
    source_local: tuple[sp.Expr, sp.Expr]
    target_local: tuple[sp.Expr, sp.Expr]
    source_inward: tuple[float, float]
    target_inward: tuple[float, float]


@dataclass(frozen=True)
class SurfaceMesh:
    name: str
    parameter_mesh: Mesh
    spatial_vertices: dict[str, tuple[sp.Expr, sp.Expr, sp.Expr]]
    vertex_classes: dict[str, str]
    edges: tuple[SurfaceEdge, ...]
    face_orientation: dict[str, int]
    vertex_normals: dict[str, np.ndarray]
    center: np.ndarray
    radius: float | None

    @property
    def vertex_ids(self):
        return tuple(sorted(set(self.vertex_classes.values())))


def load_surface_mesh(path):
    path = Path(path)
    data = json.loads(path.read_text())
    parameter = load_mesh(path.parent/data['parameter_mesh'])
    if any(not face.corners for face in parameter.faces):
        raise ValueError('Spatial faces need ordered parameter-domain corners')
    spatial = {str(key): tuple(expression(value) for value in coords)
               for key, coords in data['spatial_vertices'].items()}
    if set(spatial) != set(parameter.vertices) or any(
            len(point) != 3 or any(value.free_symbols for value in point)
            for point in spatial.values()):
        raise ValueError('Provide a constant 3D position for every parameter vertex')

    parents = {vid: vid for vid in parameter.vertices}

    def find(vid):
        while parents[vid] != vid:
            vid = parents[vid]
        return vid

    def union(first, second):
        a, b = find(first), find(second)
        if a != b:
            parents[max(a, b)] = min(a, b)

    by_edge = {edge.id: edge for edge in parameter.edges}
    face_origins = {face.id: parameter.vertices[face.local_origin]
                    for face in parameter.faces}
    face_centers = {face.id: np.mean(
        np.asarray([parameter.vertices[vid] for vid in face.corners], dtype=float), axis=0)
        for face in parameter.faces}
    physical_edges = []

    def append_edge(eid, source, target, target_parameter, pairs):
        for a, b in pairs:
            union(a, b)
            if any(sp.simplify(spatial[a][k]-spatial[b][k]) != 0 for k in range(3)):
                raise ValueError(f'Inconsistent spatial positions on edge {eid}: {a}, {b}')

        def inward(edge, fid, parameter_value):
            start = np.asarray(edge.point(parameter.vertices, parameter_value.subs(edge_t, 0)),
                               dtype=float)
            end = np.asarray(edge.point(parameter.vertices, parameter_value.subs(edge_t, 1)),
                             dtype=float)
            direction = end-start
            length = np.linalg.norm(direction)
            if length == 0:
                raise ValueError(f'Degenerate parameter edge {edge.id}')
            normal = np.array([-direction[1], direction[0]])/length
            if np.dot(normal, face_centers[fid]-start) < 0:
                normal = -normal
            return tuple(normal)

        source_parameter = edge_t
        source_point = source.point(parameter.vertices, source_parameter)
        target_point = target.point(parameter.vertices, target_parameter)
        source_origin, target_origin = face_origins[source.faces[0]], face_origins[target.faces[-1]]
        physical_edges.append(SurfaceEdge(
            eid, source.vertices, source.faces[0], target.faces[-1],
            tuple(sp.expand(source_point[k]-source_origin[k]) for k in range(2)),
            tuple(sp.expand(target_point[k]-target_origin[k]) for k in range(2)),
            inward(source, source.faces[0], source_parameter),
            inward(target, target.faces[-1], target_parameter)))

    for edge in parameter.interior_edges:
        append_edge(edge.id, edge, edge, edge_t,
                    ((vertex, vertex) for vertex in edge.vertices))
    for seam in parameter.periodic_identifications:
        source, target = by_edge[seam.source], by_edge[seam.target]
        def mapped(vertex):
            point = parameter.vertices[vertex]
            return tuple(point[k]+seam.translation[k] for k in range(2))
        pairs = []
        for vertex in source.vertices:
            hits = [other for other in target.vertices if all(
                sp.simplify(mapped(vertex)[k]-parameter.vertices[other][k]) == 0
                for k in range(2))]
            if len(hits) != 1:
                raise ValueError(f'Unmatched seam endpoint on {seam.id}')
            pairs.append((vertex, hits[0]))
        target_parameter = edge_t if pairs[0][1] == target.vertices[0] else 1-edge_t
        append_edge(seam.id, source, target, target_parameter, pairs)
    for gluing in parameter.gluings:
        source, target = by_edge[gluing.source], by_edge[gluing.target]
        pairs = []
        for vertex in source.vertices:
            mapped = gluing.matrix*sp.Matrix(parameter.vertices[vertex])+sp.Matrix(gluing.offset)
            hits = [other for other in target.vertices if all(
                sp.simplify(mapped[k]-parameter.vertices[other][k]) == 0
                for k in range(2))]
            if len(hits) != 1:
                raise ValueError(f'Unmatched gluing endpoint on {gluing.id}')
            pairs.append((vertex, hits[0]))
        target_parameter = edge_t if pairs[0][1] == target.vertices[0] else 1-edge_t
        append_edge(gluing.id, source, target, target_parameter, pairs)

    classes = {vid: find(vid) for vid in parameter.vertices}
    # Interior shared edges already have one vertex ID, while identified edges
    # each become one edge of the closed spatial mesh.
    if len(physical_edges) != len(parameter.interior_edges) + len(parameter.gluings) + len(
            parameter.periodic_identifications):
        raise ValueError('Every spatial edge must have exactly two incident faces')
    unique = {classes[vid]: np.asarray(point, dtype=float)
              for vid, point in spatial.items()}
    center = np.mean(list(unique.values()), axis=0)
    radii = np.asarray([np.linalg.norm(point-center) for point in unique.values()])
    radius = float(np.mean(radii)) if np.ptp(radii) < 1e-10 else None
    orientation = {}
    affine_normals = {}
    for face in parameter.faces:
        corner = face.corners
        q0 = np.asarray(parameter.vertices[corner[0]], dtype=float)
        p0 = np.asarray(spatial[corner[0]], dtype=float)
        vectors = []
        for vid in corner[1:]:
            q = np.asarray(parameter.vertices[vid], dtype=float)-q0
            p = np.asarray(spatial[vid], dtype=float)-p0
            vectors.append((q, p))
        for qa, pa in vectors:
            for qb, pb in vectors:
                area = np.linalg.det(np.stack((qa, qb)))
                if abs(area) > 1e-10:
                    jacobian = np.stack((pa, pb), axis=1) @ np.linalg.inv(
                        np.stack((qa, qb), axis=1))
                    normal = np.cross(jacobian[:, 0], jacobian[:, 1])
                    if radius is None:
                        orientation[face.id] = 1 if area > 0 else -1
                    else:
                        facing = np.dot(normal, np.mean(
                            [np.asarray(spatial[vid], dtype=float) for vid in corner], axis=0)-center)
                        if abs(facing) < 1e-10:
                            raise ValueError(f'Cannot orient face {face.id}')
                        orientation[face.id] = 1 if facing > 0 else -1
                    affine_normals[face.id] = orientation[face.id]*normal/np.linalg.norm(normal)
                    break
            if face.id in orientation:
                break
        if face.id not in orientation:
            raise ValueError(f'Degenerate spatial face {face.id}')
    vertex_normals = {}
    for face in parameter.faces:
        for representative in {classes[vid] for vid in face.corners}:
            vertex_normals.setdefault(representative, np.zeros(3))
            vertex_normals[representative] += affine_normals[face.id]
    for representative, normal in vertex_normals.items():
        magnitude = np.linalg.norm(normal)
        if magnitude < 1e-10:
            raise ValueError(f'No consistent reference normal at {representative}')
        vertex_normals[representative] = normal/magnitude
    return SurfaceMesh(str(data.get('name', path.stem)), parameter, spatial, classes,
                       tuple(physical_edges), orientation, vertex_normals, center, radius)
