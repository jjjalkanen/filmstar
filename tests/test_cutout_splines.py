"""Rotational seam and corner-normal checks for the polyhedron cutouts."""
from pathlib import Path
import unittest

import numpy as np
import sympy as sp

from filmstar.planar.mesh import load_mesh
from filmstar.plotting.data import load_plot_geometry
from filmstar.plotting.planar import SplineProjection
from filmstar.planar.solver import SplineSpace

ROOT = Path(__file__).resolve().parents[1]
CASES = (
    ('tetrahedron', 6, 9, 4, 3, 4),
    ('cube', 14, 19, 6, 7, 4),
    ('octahedron', 10, 17, 8, 5, 6),
    ('dodecahedron', 40, 50, 12, 20, 4),
)


class CutoutTests(unittest.TestCase):
    def test_cutout_meshes_and_c0_plots(self):
        for name, nv, ne, nf, ng, dimension in CASES:
            with self.subTest(name=name):
                path = ROOT/'examples/meshes'/f'{name}.json'
                mesh = load_mesh(path)
                self.assertEqual((len(mesh.vertices), len(mesh.edges), len(mesh.faces),
                                  len(mesh.gluings)), (nv, ne, nf, ng))
                self.assertEqual(mesh.periodic_identifications, ())
                by_id = {edge.id: edge for edge in mesh.edges}
                parents = {vertex: vertex for vertex in mesh.vertices}

                def find(vertex):
                    while parents[vertex] != vertex:
                        vertex = parents[vertex]
                    return vertex

                for gluing in mesh.gluings:
                    self.assertEqual(sp.simplify(gluing.matrix.det()), 1)
                    self.assertEqual(sp.simplify(gluing.matrix.T*gluing.matrix
                                                  - sp.eye(2)), sp.zeros(2))
                    source, target = by_id[gluing.source], by_id[gluing.target]
                    mapped = [gluing.matrix*sp.Matrix(mesh.vertices[v])
                              + sp.Matrix(gluing.offset) for v in source.vertices]
                    self.assertTrue(all(any(all(sp.simplify(point[k]
                                                         - mesh.vertices[v][k]) == 0
                                                for k in range(2))
                                            for v in target.vertices)
                                        for point in mapped))
                    if name == 'dodecahedron':
                        for vertex, point in zip(source.vertices, mapped):
                            partner = next(v for v in target.vertices if all(
                                sp.simplify(point[k]-mesh.vertices[v][k]) == 0
                                for k in range(2)))
                            parents[find(vertex)] = find(partner)
                if name == 'dodecahedron':
                    degrees = {find(vertex): set() for vertex in mesh.vertices}
                    for edge in mesh.edges:
                        a, b = (find(vertex) for vertex in edge.vertices)
                        degrees[a].add(b)
                        degrees[b].add(a)
                    self.assertEqual(len(degrees), 20)
                    self.assertEqual({len(neighbors) for neighbors in degrees.values()}, {3})
                    self.assertEqual(sum(map(len, degrees.values())) // 2, 30)
                space = SplineSpace(mesh, 1, 0)
                self.assertEqual(space.dimension, dimension)
                spline = space.spline([sp.Rational((-1)**k, k+1)
                                       for k in range(dimension)])
                self.assertTrue(spline.verify())
                projection = SplineProjection(spline, load_plot_geometry(path, mesh))
                samples = projection.surface_func(np.linspace(0, 1, 11), 0)
                heights = np.concatenate([sample[:, :, 2].ravel() for sample in samples])
                self.assertTrue(np.all(np.isfinite(heights)))
                self.assertGreater(np.ptp(heights), 0.01)

    def test_c1_normals_from_two_face_edges(self):
        for name, expected in (('tetrahedron', 4), ('cube', 4), ('octahedron', 5)):
            with self.subTest(name=name):
                mesh = load_mesh(ROOT/'examples/meshes'/f'{name}.json')
                edge_space = SplineSpace(mesh, 3, 1, vertex_normals=False)
                space = SplineSpace(mesh, 3, 1, vertex_normals=True)
                self.assertGreater(space.vertex_matrix.rows, 0)
                self.assertEqual(edge_space.dimension, expected)
                self.assertEqual(space.dimension, expected)
                spline = space.spline([sp.Rational(1, k+1) for k in range(expected)])
                self.assertTrue(spline.verify())
                self.assertEqual(space.vertex_matrix*space.kernel,
                                 sp.zeros(space.vertex_matrix.rows, expected))

    def test_dodecahedron_c1_linear_space(self):
        mesh = load_mesh(ROOT/'examples/meshes/dodecahedron.json')
        space = SplineSpace(mesh, 1, 1)
        self.assertGreater(space.vertex_matrix.rows, 0)
        self.assertEqual(space.dimension, 1)
        self.assertEqual(space.vertex_matrix*space.kernel,
                         sp.zeros(space.vertex_matrix.rows, 1))
        self.assertTrue(space.spline([sp.Rational(1, 2)]).verify())


if __name__ == '__main__':
    unittest.main()
