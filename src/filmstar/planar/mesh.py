"""Mesh topology and edge equations for the exact spline calculation.

Plotting grids in the JSON file are intentionally ignored here. A spline
space depends on faces, their adjacent edges, and the edge polynomials.
"""
from dataclasses import dataclass
import json
from pathlib import Path

import sympy as sp

x, y, edge_t = sp.symbols('x y t')
_ALLOWED = {'x': x, 'y': y, 't': edge_t,
            'sqrt': sp.sqrt, 'Abs': sp.Abs, 'Rational': sp.Rational}


def expression(value):
    result = sp.sympify(value, locals=_ALLOWED)
    if result.free_symbols - {x, y, edge_t}:
        raise ValueError(f'Unexpected symbol in mesh expression: {value}')
    return result


@dataclass(frozen=True)
class Edge:
    id: str
    vertices: tuple[str, str]
    polynomial: sp.Expr
    faces: tuple[str, ...]
    curve: tuple[sp.Expr, sp.Expr]

    def point(self, vertices, parameter):
        """Physical point on this edge at parameter 0 <= t <= 1."""
        start = vertices[self.vertices[0]]
        return tuple(start[k] + self.curve[k].subs(edge_t, parameter)
                     for k in range(2))


@dataclass(frozen=True)
class Face:
    id: str
    local_origin: str
    corners: tuple[str, ...]


@dataclass(frozen=True)
class PeriodicIdentification:
    """Map a source boundary edge to a target edge by a translation."""
    id: str
    source: str
    target: str
    translation: tuple[sp.Expr, sp.Expr]


@dataclass(frozen=True)
class Gluing:
    """An orientation-preserving rigid map from source to target coordinates."""
    id: str
    source: str
    target: str
    matrix: sp.ImmutableMatrix
    offset: tuple[sp.Expr, sp.Expr]


@dataclass(frozen=True)
class Mesh:
    name: str
    vertices: dict[str, tuple[sp.Expr, sp.Expr]]
    edges: tuple[Edge, ...]
    faces: tuple[Face, ...]
    periodic_identifications: tuple[PeriodicIdentification, ...]
    gluings: tuple[Gluing, ...]

    @property
    def interior_edges(self):
        return tuple(edge for edge in self.edges if len(edge.faces) == 2)


def load_mesh(path):
    """Load and validate topology; no face grid is needed or read."""
    data = json.loads(Path(path).read_text())
    vertices = {str(k): tuple(expression(a) for a in xy)
                for k, xy in data['vertices'].items()}
    if any(len(point) != 2 or any(value.free_symbols for value in point)
           for point in vertices.values()):
        raise ValueError('Vertices must have two constant coordinates')
    faces = tuple(Face(str(f['id']), str(f['local_origin']),
                       tuple(map(str, f.get('corners', ())))) for f in data['faces'])
    face_ids = [face.id for face in faces]
    if not faces or len(face_ids) != len(set(face_ids)):
        raise ValueError('Faces must have unique IDs')
    edges = tuple(Edge(str(edge['id']), tuple(map(str, edge['vertices'])),
                       expression(edge['polynomial']), tuple(map(str, edge['faces'])),
                       (expression(edge['curve']['dx']), expression(edge['curve']['dy'])))
                  for edge in data['edges'])
    if len({edge.id for edge in edges}) != len(edges):
        raise ValueError('Edges must have unique IDs')
    for edge in edges:
        if (len(edge.vertices) != 2 or len(set(edge.vertices)) != 2
                or not 1 <= len(edge.faces) <= 2 or len(set(edge.faces)) != len(edge.faces)
                or any(fid not in face_ids for fid in edge.faces)
                or any(vertex not in vertices for vertex in edge.vertices)):
            raise ValueError(f'Invalid incidence for edge {edge.id}')
        poly = sp.Poly(edge.polynomial, x, y)
        if poly.is_zero or poly.total_degree() < 1 or edge.polynomial.has(edge_t):
            raise ValueError(f'Invalid defining polynomial for edge {edge.id}')
        for vertex_id in edge.vertices:
            vx, vy = vertices[vertex_id]
            if sp.simplify(poly.as_expr().subs({x: vx, y: vy})) != 0:
                raise ValueError(f'Edge {edge.id} misses vertex {vertex_id}')
        if any(component.free_symbols - {edge_t} for component in edge.curve):
            raise ValueError(f'Edge {edge.id} curve must depend only on t')
        for parameter, vertex_id in ((0, edge.vertices[0]), (1, edge.vertices[1])):
            if any(sp.simplify(a-b) != 0 for a, b in zip(edge.point(vertices, parameter),
                                                         vertices[vertex_id])):
                raise ValueError(f'Edge {edge.id} curve misses endpoint {vertex_id}')
        px, py = edge.point(vertices, edge_t)
        if sp.simplify(edge.polynomial.subs({x: px, y: py}, simultaneous=True)) != 0:
            raise ValueError(f'Edge {edge.id} curve misses its defining equation')
    for face in faces:
        if sum(face.id in edge.faces for edge in edges) < 3:
            raise ValueError(f'Face {face.id} has fewer than three incident edges')
        incident_vertices = {v for edge in edges if face.id in edge.faces
                             for v in edge.vertices}
        if face.local_origin not in incident_vertices:
            raise ValueError(f'Face {face.id} local origin is not a corner')
        if face.corners:
            if len(face.corners) < 3 or len(set(face.corners)) != len(face.corners):
                raise ValueError(f'Face {face.id} has invalid ordered corners')
            if set(face.corners) != incident_vertices:
                raise ValueError(f'Face {face.id} corners disagree with its edges')
            actual_edges = {frozenset(edge.vertices) for edge in edges if face.id in edge.faces}
            cycle_edges = {frozenset((face.corners[k], face.corners[(k+1)%len(face.corners)]))
                           for k in range(len(face.corners))}
            if actual_edges != cycle_edges:
                raise ValueError(f'Face {face.id} corners do not follow its edges')
    identifications = tuple(
        PeriodicIdentification(str(item['id']), str(item['source']),
                               str(item['target']),
                               tuple(expression(value) for value in item['translation']))
        for item in data.get('periodic_identifications', ()))
    if identifications:
        ids = [item.id for item in identifications]
        if len(ids) != len(set(ids)) or set(ids) & {edge.id for edge in edges}:
            raise ValueError('Periodic identification IDs must be unique and distinct from edge IDs')
        by_id = {edge.id: edge for edge in edges}
        paired_edges = []
        for item in identifications:
            if item.source not in by_id or item.target not in by_id:
                raise ValueError(f'Unknown edge in periodic identification {item.id}')
            source, target = by_id[item.source], by_id[item.target]
            if len(source.faces) != 1 or len(target.faces) != 1 or source.id == target.id:
                raise ValueError(f'Periodic identification {item.id} must pair distinct boundary edges')
            if len(item.translation) != 2 or any(value.free_symbols for value in item.translation):
                raise ValueError(f'Periodic identification {item.id} needs a constant 2D translation')
            dx, dy = item.translation
            shifted = {(sp.simplify(vertices[v][0]+dx), sp.simplify(vertices[v][1]+dy))
                       for v in source.vertices}
            destination = {vertices[v] for v in target.vertices}
            if shifted != destination:
                raise ValueError(f'Periodic identification {item.id} does not map edge endpoints')
            sx, sy = source.point(vertices, edge_t)
            target_parameter = edge_t if all(
                sp.simplify(a+b-c) == 0 for a, b, c in
                zip(source.point(vertices, 0), item.translation,
                    target.point(vertices, 0))) else 1-edge_t
            tx, ty = target.point(vertices, target_parameter)
            if (sp.simplify(sx+dx-tx) != 0 or sp.simplify(sy+dy-ty) != 0):
                raise ValueError(f'Periodic identification {item.id} does not map edge parameterizations')
            translated_curve = sp.Poly(target.polynomial.subs({x: x+dx, y: y+dy}, simultaneous=True), x, y)
            source_curve = sp.Poly(source.polynomial, x, y)
            if (sp.expand(translated_curve.as_expr()*source_curve.LC()
                          - source_curve.as_expr()*translated_curve.LC()) != 0):
                raise ValueError(f'Periodic identification {item.id} does not map edge curves')
            paired_edges.extend((source.id, target.id))
        boundary = {edge.id for edge in edges if len(edge.faces) == 1}
        if len(paired_edges) != len(set(paired_edges)) or set(paired_edges) != boundary:
            raise ValueError('A periodic mesh must pair every boundary edge exactly once')
        translations = sp.Matrix([list(item.translation) for item in identifications])
        if translations.rank() != 2:
            raise ValueError('A torus needs two independent periods')
    gluings = tuple(_load_gluing(item, vertices, edges, faces)
                    for item in data.get('edge_identifications', ()))
    if gluings:
        by_id = {edge.id: edge for edge in edges}
        paired = [eid for item in gluings for eid in (item.source, item.target)]
        boundary = {edge.id for edge in edges if len(edge.faces) == 1}
        if (len(paired) != len(set(paired)) or set(paired) != boundary
                or len({item.id for item in gluings}) != len(gluings)
                or {item.id for item in gluings} & {edge.id for edge in edges}):
            raise ValueError('Gluings must pair every boundary edge exactly once')
        if identifications:
            raise ValueError('Use either rigid edge identifications or periodic translations')
        for item in gluings:
            source, target = by_id[item.source], by_id[item.target]
            if len(source.faces) != 1 or len(target.faces) != 1:
                raise ValueError(f'Gluing {item.id} must pair boundary edges')
    return Mesh(str(data.get('name', Path(path).stem)), vertices, edges, faces,
                identifications, gluings)


def _load_gluing(item, vertices, edges, faces):
    """Derive the direct isometry from oppositely oriented boundary edges."""
    by_id = {edge.id: edge for edge in edges}
    source_id, target_id = str(item['source']), str(item['target'])
    if source_id not in by_id or target_id not in by_id:
        raise ValueError('Gluing refers to an unknown edge')
    source, target = by_id[source_id], by_id[target_id]
    if source_id == target_id or len(source.faces) != 1 or len(target.faces) != 1:
        raise ValueError('Gluing needs two distinct boundary edges')

    def oriented(edge):
        face = next(f for f in faces if f.id == edge.faces[0])
        if not face.corners:
            raise ValueError(f'Glued face {face.id} needs ordered corners')
        cycle = list(face.corners)
        area = sp.expand(sum(
            vertices[cycle[k]][0]*vertices[cycle[(k+1)%len(cycle)]][1]
            - vertices[cycle[(k+1)%len(cycle)]][0]*vertices[cycle[k]][1]
            for k in range(len(cycle))))
        if float(area) < 0:
            cycle.reverse()
        for k in range(len(cycle)):
            pair = (cycle[k], cycle[(k+1)%len(cycle)])
            if set(pair) == set(edge.vertices):
                return pair
        raise ValueError(f'Edge {edge.id} missing from face cycle')

    sa, sb = oriented(source)
    ta, tb = oriented(target)
    p, q = sp.Matrix(vertices[sa]), sp.Matrix(vertices[sb])
    target_start, target_end = sp.Matrix(vertices[tb]), sp.Matrix(vertices[ta])
    a, b = q-p, target_end-target_start
    length_squared = sp.expand(a.dot(a))
    if sp.simplify(length_squared-b.dot(b)) != 0:
        raise ValueError(f'Gluing {item["id"]} pairs edges of unequal length')
    cosine = sp.simplify(a.dot(b)/length_squared)
    sine = sp.simplify((a[0]*b[1]-a[1]*b[0])/length_squared)
    rotation = sp.ImmutableMatrix([[cosine, -sine], [sine, cosine]])
    offset = tuple(sp.simplify(value) for value in target_start-rotation*p)
    return Gluing(str(item['id']), source_id, target_id, rotation, offset)
