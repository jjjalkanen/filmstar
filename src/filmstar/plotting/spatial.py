"""Sample solved spatial patches for the shared 3D plotter."""
import numpy as np

from filmstar.plotting.data import PolygonGrid
from filmstar.spatial.solver import SurfaceSolution


class SurfaceProjection:
    def __init__(self, solution: SurfaceSolution):
        if not isinstance(solution, SurfaceSolution):
            raise TypeError('Expected a solved spatial surface')
        self.solution = solution
        mesh = solution.space.mesh.parameter_mesh
        self.grids = tuple(PolygonGrid(face.id,
                                       tuple(mesh.vertices[vid] for vid in face.corners))
                           for face in mesh.faces)

    @property
    def title(self):
        space = self.solution.space
        return rf'{space.mesh.name}: $C^{{{space.smoothness}}}_{{{space.degree}}}$ spatial surface'

    def surface_func(self, tau, alpha):
        c, sine = np.cos(alpha), np.sin(alpha)
        surfaces = []
        for grid in self.grids:
            xx, yy = grid.sample(tau)
            xyz = self.solution.evaluate_face(grid.face_id, xx, yy)
            rotated = np.empty_like(xyz)
            rotated[:, :, 0] = c*xyz[:, :, 0]-sine*xyz[:, :, 1]
            rotated[:, :, 1] = sine*xyz[:, :, 0]+c*xyz[:, :, 1]
            rotated[:, :, 2] = xyz[:, :, 2]
            surfaces.append(rotated)
        return surfaces
