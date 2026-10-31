"""Generate 3D plotting data from one solved spline and its face grids."""
import numpy as np

from filmstar.plotting.data import PlotGeometry
from filmstar.planar.solver import Spline


class SplineProjection:
    def __init__(self, spline: Spline, geometry: PlotGeometry):
        if not isinstance(spline, Spline) or not isinstance(geometry, PlotGeometry):
            raise TypeError('Expected a solved Spline and PlotGeometry')
        if tuple(spline.pieces) != geometry.face_ids:
            raise ValueError('Plot geometry does not match spline faces')
        self.spline, self.geometry = spline, geometry

    @property
    def title(self):
        space = self.spline.space
        return rf'{self.geometry.name}: $C^{{{space.smoothness}}}_{{{space.degree}}}$ spline'

    def surface_func(self, tau, alpha):
        c, s = np.cos(alpha), np.sin(alpha)
        surfaces = []
        for grid in self.geometry.grids:
            xx, yy = grid.sample(tau)
            zz = self.spline.evaluate_face(grid.face_id, xx, yy)
            surfaces.append(np.stack((c*xx-s*yy, s*xx+c*yy, zz), axis=-1))
        return surfaces
