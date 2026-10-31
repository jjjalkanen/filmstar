"""Topology-only solving and separate plotting data across three meshes."""
import json
from dataclasses import replace
import subprocess
import tempfile
from pathlib import Path
import unittest

import numpy as np
import sympy as sp

from filmstar.planar.mesh import load_mesh, x, y
from filmstar.plotting.data import load_plot_geometry
from filmstar.plotting.planar import SplineProjection
from filmstar.planar.solver import SplineSpace, s, t
from filmstar.cli.splines import DEMO_COEFFICIENTS

ROOT = Path(__file__).resolve().parents[1]
CASES = (
    ('paper', 8, 9, 6, 26, (3, 0)),
    ('triangles', 2, 1, 4, 16, (3, 0)),
    ('squares', 5, 8, 4, 23, (3, 0)),
)
PERIODIC_CASES = (
    ('paper_torus', 3, 0, 3, 18),
    ('paper_torus', 5, 1, 3, 12),
    ('triangles_torus', 3, 0, 2, 9),
    ('squares_torus', 3, 0, 2, 18),
)


def example_spline(space):
    return space.spline([DEMO_COEFFICIENTS[k % len(DEMO_COEFFICIENTS)]
                         for k in range(space.dimension)])


class SplineTests(unittest.TestCase):
    def test_all_mesh_files(self):
        for name, faces, interior, boundary, dimension, (degree, smoothness) in CASES:
            with self.subTest(mesh=name):
                path = ROOT/'examples/meshes'/f'{name}.json'
                mesh = load_mesh(path)
                self.assertEqual(len(mesh.faces), faces)
                self.assertEqual(len(mesh.interior_edges), interior)
                self.assertEqual(len(mesh.edges)-interior, boundary)
                self.assertEqual(mesh.periodic_identifications, ())
                for edge in mesh.edges:
                    self.assertEqual(edge.point(mesh.vertices, 0),
                                     mesh.vertices[edge.vertices[0]])
                    self.assertEqual(edge.point(mesh.vertices, 1),
                                     mesh.vertices[edge.vertices[1]])
                    px, py = edge.point(mesh.vertices, sp.Rational(1, 2))
                    self.assertEqual(sp.simplify(edge.polynomial.subs(
                        {x: px, y: py}, simultaneous=True)), 0)
                space = SplineSpace(mesh, degree, smoothness)
                self.assertEqual(space.dimension, dimension)
                self.assertTrue(all(sp.simplify(value) == 0
                                    for value in space.matrix*space.kernel))
                spline = example_spline(space)
                self.assertTrue(spline.verify())
                self.assertGreater(len(set(spline.pieces.values())), 1)
                for face in mesh.faces:
                    ox, oy = mesh.vertices[face.local_origin]
                    self.assertEqual(sp.simplify(spline.pieces[face.id].subs({s: 0, t: 0})
                                                 - spline.physical_polynomial(face.id).subs(
                                                     {x: ox, y: oy}, simultaneous=True)), 0)
                geometry = load_plot_geometry(path, mesh)
                projection = SplineProjection(spline, geometry)
                surfaces = projection.surface_func(np.linspace(0, 1, 21), 0)
                self.assertEqual(len(surfaces), faces)
                points = np.concatenate([surface.reshape(-1, 3) for surface in surfaces])
                self.assertTrue(np.all(np.isfinite(points)))
                self.assertGreater(np.ptp(points[:, 2]), 0.1)
                with tempfile.TemporaryDirectory() as tmp:
                    output = Path(tmp)/'surface.png'
                    completed = subprocess.run(
                        [str(ROOT/'run.sh'), '--mesh', str(path), '--deg', str(degree),
                         '--smooth', str(smoothness), '--save', str(output)],
                        cwd=ROOT, capture_output=True, text=True, check=True)
                    self.assertIn(f'dimension {dimension}', completed.stdout)
                    self.assertGreater(output.stat().st_size, 10000)

    def test_periodic_mesh_files(self):
        for name, degree, smoothness, seam_count, dimension in PERIODIC_CASES:
            with self.subTest(mesh=name, degree=degree, smoothness=smoothness):
                path = ROOT/'examples/meshes'/f'{name}.json'
                mesh = load_mesh(path)
                self.assertEqual(len(mesh.periodic_identifications), seam_count)
                space = SplineSpace(mesh, degree, smoothness)
                self.assertEqual(space.dimension, dimension)
                self.assertTrue(all(sp.simplify(value) == 0
                                    for value in space.matrix*space.kernel))
                spline = example_spline(space)
                self.assertTrue(spline.verify())
                self.assertGreater(len(set(spline.pieces.values())), 1)
                by_id = {edge.id: edge for edge in mesh.edges}
                for seam in mesh.periodic_identifications:
                    source = by_id[seam.source]
                    target = by_id[seam.target]
                    dx, dy = seam.translation
                    # Check traces independently from the stored cofactors.
                    delta = sp.expand(spline.physical_polynomial(target.faces[0]).subs(
                        {x: x+dx, y: y+dy}, simultaneous=True)
                        - spline.physical_polynomial(source.faces[0]))
                    for total in range(smoothness+1):
                        for nx in range(total+1):
                            derivative = sp.diff(delta, x, nx, y, total-nx)
                            if sp.Poly(source.polynomial, x, y).degree(y) == 1:
                                trace = derivative.subs(y, sp.solve(source.polynomial, y)[0])
                            else:
                                trace = derivative.subs(x, sp.solve(source.polynomial, x)[0])
                            self.assertEqual(sp.expand(trace), 0)
                geometry = load_plot_geometry(path, mesh)
                surfaces = SplineProjection(spline, geometry).surface_func(
                    np.linspace(0, 1, 13), 0)
                self.assertEqual(len(surfaces), len(mesh.faces))
                self.assertTrue(all(np.all(np.isfinite(surface)) for surface in surfaces))
                with tempfile.TemporaryDirectory() as tmp:
                    output = Path(tmp)/'surface.png'
                    completed = subprocess.run(
                        [str(ROOT/'run.sh'), '--mesh', str(path), '--deg', str(degree),
                         '--smooth', str(smoothness), '--save', str(output)],
                        cwd=ROOT, capture_output=True, text=True, check=True)
                    self.assertIn(f'dimension {dimension}', completed.stdout)
                    self.assertGreater(output.stat().st_size, 10000)

    def test_solver_ignores_plotting_grids(self):
        source = json.loads((ROOT/'examples/meshes/triangles.json').read_text())
        for face in source['faces']:
            face.pop('grid')
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'topology_only.json'
            path.write_text(json.dumps(source))
            mesh = load_mesh(path)
            space = SplineSpace(mesh, 3, 0)
            self.assertEqual(space.dimension, 16)
            self.assertTrue(example_spline(space).verify())
            with self.assertRaises(KeyError):
                load_plot_geometry(path, mesh)

    def test_paper_dimensions_and_all_basis_traces(self):
        mesh = load_mesh(ROOT/'examples/meshes/paper.json')
        for degree, smoothness, shape, rank, dimension in (
                (3, 0, (90, 110), 84, 26), (5, 1, (189, 202), 168, 34)):
            with self.subTest(degree=degree, smoothness=smoothness):
                space = SplineSpace(mesh, degree, smoothness)
                self.assertEqual(space.matrix.shape, shape)
                self.assertEqual(space.rank, rank)
                self.assertEqual(space.dimension, dimension)
                self.assertEqual(space.kernel[:8*len(space.monomials), :].rank(), dimension)
                self.assertTrue(example_spline(space).verify())
                for pieces in space.basis:
                    for edge in mesh.interior_edges:
                        physical = {}
                        for face in mesh.faces:
                            ox, oy = mesh.vertices[face.local_origin]
                            physical[face.id] = pieces[face.id].subs(
                                {s: x-ox, t: y-oy}, simultaneous=True)
                        delta = physical[edge.faces[0]] - physical[edge.faces[1]]
                        edge_y = sp.solve(edge.polynomial, y)[0]
                        for total in range(smoothness+1):
                            for dx in range(total+1):
                                trace = sp.diff(delta, x, dx, y, total-dx).subs(y, edge_y)
                                self.assertEqual(sp.expand(trace), 0)

    def test_invalid_edge_endpoint_is_rejected(self):
        source = json.loads((ROOT/'examples/meshes/triangles.json').read_text())
        source['edges'][0]['polynomial'] = 'y-1'
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'bad.json'
            path.write_text(json.dumps(source))
            with self.assertRaisesRegex(ValueError, 'misses vertex'):
                load_mesh(path)

    def test_invalid_periodic_pairing_is_rejected(self):
        source = json.loads((ROOT/'examples/meshes/triangles_torus.json').read_text())
        source['periodic_identifications'][0]['translation'] = ['0', '-sqrt(3)']
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'bad_period.json'
            path.write_text(json.dumps(source))
            with self.assertRaisesRegex(ValueError, 'does not map edge endpoints'):
                load_mesh(path)

    def test_large_global_translation_preserves_local_system(self):
        mesh = load_mesh(ROOT/'examples/meshes/paper_torus.json')
        dx, dy = sp.Integer(10**9), sp.Integer(-10**9)
        translated = replace(
            mesh,
            vertices={key: (point[0]+dx, point[1]+dy)
                      for key, point in mesh.vertices.items()},
            edges=tuple(replace(edge, polynomial=sp.expand(edge.polynomial.subs(
                {x: x-dx, y: y-dy}, simultaneous=True))) for edge in mesh.edges))
        original_space = SplineSpace(mesh, 3, 0)
        translated_space = SplineSpace(translated, 3, 0)
        self.assertEqual(original_space.matrix, translated_space.matrix)
        self.assertEqual(original_space.kernel, translated_space.kernel)
        original = example_spline(original_space)
        shifted = example_spline(translated_space)
        self.assertEqual(original.pieces, shifted.pieces)
        for face in mesh.faces:
            xx = np.array([0.123, -0.321])
            yy = np.array([0.234, -0.432])
            np.testing.assert_allclose(
                shifted.evaluate_face(face.id, xx+float(dx), yy+float(dy)),
                original.evaluate_face(face.id, xx, yy), atol=1e-6, rtol=1e-6)

    def test_invalid_local_geometry_is_rejected(self):
        source = json.loads((ROOT/'examples/meshes/triangles.json').read_text())
        source['edges'][0]['curve']['dy'] = 't'
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'bad_curve.json'
            path.write_text(json.dumps(source))
            with self.assertRaisesRegex(ValueError, 'curve misses endpoint'):
                load_mesh(path)


if __name__ == '__main__':
    unittest.main()
