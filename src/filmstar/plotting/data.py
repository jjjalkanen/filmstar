"""Face sampling grids loaded separately from the spline mesh topology."""
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
import sympy as sp

from filmstar.planar.mesh import Mesh, expression, y


@dataclass(frozen=True)
class FaceGrid:
    face_id: str
    y_range: tuple[sp.Expr, sp.Expr]
    x_range: tuple[sp.Expr, sp.Expr]

    def sample(self, tau):
        tau = np.asarray(tau, dtype=float)
        if (tau.ndim != 1 or len(tau) < 2 or not np.all(np.isfinite(tau))
                or np.any((tau < 0) | (tau > 1))):
            raise ValueError('tau must be a finite 1D array in [0,1] of length >= 2')
        u, v = np.meshgrid(tau, tau)
        ymin, ymax = (float(z) for z in self.y_range)
        yy = ymin + (ymax-ymin)*v
        lo = np.asarray(sp.lambdify(y, self.x_range[0], 'numpy')(yy), dtype=float)
        hi = np.asarray(sp.lambdify(y, self.x_range[1], 'numpy')(yy), dtype=float)
        return np.broadcast_to(lo, yy.shape) + u*(hi-lo), yy


@dataclass(frozen=True)
class PolygonGrid:
    """Sample a convex polygon by horizontal slices between its two sides."""
    face_id: str
    corners: tuple[tuple[sp.Expr, sp.Expr], ...]

    def sample(self, tau):
        tau = np.asarray(tau, dtype=float)
        if (tau.ndim != 1 or len(tau) < 2 or not np.all(np.isfinite(tau))
                or np.any((tau < 0) | (tau > 1))):
            raise ValueError('tau must be a finite 1D array in [0,1] of length >= 2')
        points = np.asarray(self.corners, dtype=float)
        ymin, ymax = points[:, 1].min(), points[:, 1].max()
        levels = ymin+(ymax-ymin)*tau
        bounds = []
        for level in levels:
            intersections = []
            for a, b in zip(points, np.roll(points, -1, axis=0)):
                if abs(a[1]-b[1]) < 1e-12:
                    if abs(level-a[1]) < 1e-10:
                        intersections.extend((a[0], b[0]))
                elif min(a[1], b[1])-1e-10 <= level <= max(a[1], b[1])+1e-10:
                    intersections.append(a[0]+(level-a[1])*(b[0]-a[0])/(b[1]-a[1]))
            if not intersections:
                raise ValueError(f'No polygon slice for {self.face_id} at y={level}')
            bounds.append((min(intersections), max(intersections)))
        bounds = np.asarray(bounds)
        xx = bounds[:, 0, None]+(bounds[:, 1]-bounds[:, 0])[:, None]*tau[None, :]
        yy = np.broadcast_to(levels[:, None], xx.shape)
        return xx, yy


@dataclass(frozen=True)
class PlotGeometry:
    name: str
    grids: tuple[FaceGrid | PolygonGrid, ...]

    @property
    def face_ids(self):
        return tuple(grid.face_id for grid in self.grids)


def load_plot_geometry(path, mesh: Mesh):
    """Load face grids and check that they describe the supplied solved mesh."""
    if not isinstance(mesh, Mesh):
        raise TypeError('mesh must be a loaded Mesh')
    data = json.loads(Path(path).read_text())
    by_id = {face.id: face for face in mesh.faces}
    grids = []
    for face_data in data['faces']:
        fid = str(face_data['id'])
        if 'grid' in face_data:
            grid = FaceGrid(fid,
                            tuple(expression(value) for value in face_data['grid']['y']),
                            tuple(expression(value) for value in face_data['grid']['x']))
        elif by_id[fid].corners:
            grid = PolygonGrid(fid, tuple(mesh.vertices[vertex]
                                          for vertex in by_id[fid].corners))
        else:
            raise KeyError(f'Face {fid} needs a grid or ordered corners')
        grids.append(grid)
    grids = tuple(grids)
    if tuple(grid.face_id for grid in grids) != tuple(face.id for face in mesh.faces):
        raise ValueError('Plot grids and spline mesh have different faces or face order')
    for grid in grids:
        if isinstance(grid, PolygonGrid):
            grid.sample(np.linspace(0, 1, 13))
            continue
        if len(grid.y_range) != 2 or len(grid.x_range) != 2:
            raise ValueError(f'Invalid ranges for face {grid.face_id}')
        if any(value.free_symbols for value in grid.y_range):
            raise ValueError(f'y range must be constant for face {grid.face_id}')
        if any(value.free_symbols - {y} for value in grid.x_range):
            raise ValueError(f'x range must depend only on y for face {grid.face_id}')
        if float(grid.y_range[0]) >= float(grid.y_range[1]):
            raise ValueError(f'Reversed y range for face {grid.face_id}')
        face = next(face for face in mesh.faces if face.id == grid.face_id)
        ox, oy = mesh.vertices[face.local_origin]
        if (sp.simplify(grid.x_range[0].subs(y, grid.y_range[0])-ox) != 0
                or sp.simplify(grid.y_range[0]-oy) != 0):
            raise ValueError(f'Face {grid.face_id} grid origin differs from its local corner')
        xx, yy = grid.sample(np.linspace(0, 1, 13))
        if not np.all(np.isfinite(xx)) or not np.all(np.isfinite(yy)):
            raise ValueError(f'Nonfinite grid for face {grid.face_id}')
        if np.any(np.diff(xx, axis=1) < -1e-10):
            raise ValueError(f'Reversed x range for face {grid.face_id}')
    return PlotGeometry(mesh.name, grids)
