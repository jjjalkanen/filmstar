"""Spatial polynomial patches over four closed polyhedral meshes."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np

from filmstar.spatial.mesh import load_surface_mesh
from filmstar.plotting.spatial import SurfaceProjection
from filmstar.spatial.solver import SurfaceSolution, SurfaceSpace


ROOT = Path(__file__).resolve().parents[1]
CASES = (('tetrahedron', 4, 6, 4, 5),
         ('cube', 8, 12, 6, 5),
         ('octahedron', 6, 12, 8, 5),
         ('dodecahedron', 20, 30, 12, 6))


class SpatialSurfaceTests(unittest.TestCase):
    def test_closed_meshes_and_linear_c0(self):
        for name, vertices, edges, faces, _degree in CASES:
            with self.subTest(name=name):
                mesh = load_surface_mesh(ROOT/'examples/surface_meshes'/f'{name}.json')
                self.assertEqual((len(mesh.vertex_ids), len(mesh.edges),
                                  len(mesh.parameter_mesh.faces)),
                                 (vertices, edges, faces))
                self.assertEqual(vertices-edges+faces, 2)
                space = SurfaceSpace(mesh, 1, 0)
                surface = space.c0_surface()
                self.assertLess(np.max(np.abs(space.c0_matrix@surface.coefficients.T
                                              -space.c0_targets)), 1e-9)
                self.assertEqual(len(SurfaceProjection(surface).surface_func(
                    np.linspace(0, 1, 5), 0)), faces)
                if name == 'dodecahedron':
                    degrees = {vertex: set() for vertex in mesh.vertex_ids}
                    for edge in mesh.edges:
                        a, b = (mesh.vertex_classes[v] for v in edge.source_vertices)
                        degrees[a].add(b)
                        degrees[b].add(a)
                    self.assertEqual({len(neighbors) for neighbors in degrees.values()}, {3})

    def test_saved_c1_surfaces(self):
        for name, _vertices, _edges, faces, degree in CASES:
            with self.subTest(name=name):
                mesh = load_surface_mesh(ROOT/'examples/surface_meshes'/f'{name}.json')
                space = SurfaceSpace(mesh, degree, 1)
                saved = json.loads((ROOT/'examples/surface_meshes'/f'{name}_c1_solution.json').read_text())
                self.assertEqual((saved['degree'], saved['smoothness']), (degree, 1))
                coefficients = np.asarray(saved['coefficients'], dtype=float)
                self.assertEqual(coefficients.shape, (3, space.columns))
                solution = SurfaceSolution(space, coefficients,
                                           (coefficients-space.particular)@space.nullspace)
                self.assertLess(np.max(np.abs(space.c0_matrix@coefficients.T
                                              -space.c0_targets)), 1e-8)
                edge_error, normal_error, minimum = space.c1_errors(coefficients)
                self.assertLess(edge_error, 1e-6)
                self.assertLess(normal_error, 1e-5)
                self.assertGreater(minimum, 1e-8)
                self.assertLess(solution.independent_edge_error(), 1e-6)
                metrics = solution.geometry_metrics()
                self.assertGreater(metrics['minimum_normal_ratio'], 1e-3)
                self.assertGreater(metrics['minimum_outward_cosine'], 0.05)
                self.assertLess(metrics['maximum_radius_deviation'], 0.3)
                surfaces = SurfaceProjection(solution).surface_func(np.linspace(0, 1, 9), 0)
                self.assertEqual(len(surfaces), faces)
                self.assertTrue(all(np.all(np.isfinite(values)) for values in surfaces))

    def test_fresh_c1_solve_and_command_exit(self):
        mesh = load_surface_mesh(ROOT/'examples/surface_meshes/tetrahedron.json')
        space = SurfaceSpace(mesh, 5, 1)
        solution = space.solve_c1(max_nfev=600)
        self.assertLess(solution.independent_edge_error(), 1e-6)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/'surface.png'
            command = [str(ROOT/'run_surface.sh'), '--mesh',
                       str(ROOT/'examples/surface_meshes/tetrahedron.json'), '--deg', '5',
                       '--smooth', '1', '--solution',
                       str(ROOT/'examples/surface_meshes/tetrahedron_c1_solution.json'),
                       '--save', str(output)]
            completed = subprocess.run(command, cwd=ROOT, capture_output=True,
                                       text=True, timeout=30, check=True)
            self.assertIn('Saved', completed.stdout)
            self.assertGreater(output.stat().st_size, 10000)

    def test_torus_c1_without_spherical_reference(self):
        mesh = load_surface_mesh(ROOT/'examples/surface_meshes/torus.json')
        self.assertEqual((len(mesh.vertex_ids), len(mesh.edges),
                          len(mesh.parameter_mesh.faces)), (9, 18, 9))
        self.assertIsNone(mesh.radius)
        space = SurfaceSpace(mesh, 6, 1)
        saved = json.loads((ROOT/'examples/surface_meshes/torus_c1_solution.json').read_text())
        coefficients = np.asarray(saved['coefficients'])
        surface = SurfaceSolution(space, coefficients,
                                  (coefficients-space.particular)@space.nullspace)
        self.assertLess(np.max(np.abs(space.c0_matrix@coefficients.T
                                      -space.c0_targets)), 1e-9)
        self.assertLess(space.c1_errors(coefficients)[0], 1e-6)
        self.assertLess(space.c1_errors(coefficients)[1], 1e-5)
        self.assertLess(surface.independent_edge_error(), 1e-6)
        self.assertGreater(surface.geometry_metrics()['minimum_normal_ratio'], 1e-3)
        self.assertEqual(len(SurfaceProjection(surface).surface_func(
            np.linspace(0, 1, 9), 0)), 9)


if __name__ == '__main__':
    unittest.main()
