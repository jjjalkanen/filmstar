"""Solve and plot a spline from a JSON mesh, interactively or to a PNG."""
import argparse
from pathlib import Path

import numpy as np
import sympy as sp

from filmstar.planar.mesh import load_mesh
from filmstar.plotting.data import load_plot_geometry
from filmstar.plotting.planar import SplineProjection
from filmstar.planar.solver import SplineSpace

# Explicit, reproducible example coefficients. SplineSpace itself selects none.
DEMO_COEFFICIENTS = tuple(map(sp.Rational, ('-1', '-2/3', '-1/2', '-1/3',
                                            '1/3', '1/2', '2/3', '1')))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=Path, required=True, help='JSON mesh file')
    parser.add_argument('--deg', type=int, required=True, help='maximum polynomial degree')
    parser.add_argument('--smooth', type=int, required=True, help='continuity order')
    parser.add_argument('--save', type=Path, metavar='PNG', help='write a PNG instead of opening a window')
    args = parser.parse_args()

    import matplotlib
    matplotlib.use('Agg' if args.save else 'qtagg')

    mesh = load_mesh(args.mesh)
    space = SplineSpace(mesh, args.deg, args.smooth,
                        vertex_normals=bool(mesh.gluings and args.smooth >= 1))
    coefficients = [DEMO_COEFFICIENTS[k % len(DEMO_COEFFICIENTS)]
                    for k in range(space.dimension)]
    spline = space.spline(coefficients)
    if not spline.verify():
        raise RuntimeError('Polynomial matching identities failed')
    geometry = load_plot_geometry(args.mesh, mesh)
    projection = SplineProjection(spline, geometry)
    print(f'{mesh.name}: {len(mesh.faces)} faces, {len(mesh.interior_edges)} interior edges, '
          f'{len(mesh.periodic_identifications)} periodic identifications, '
          f'{len(mesh.gluings)} rigid edge identifications, '
          f'C^{args.smooth}_{args.deg} dimension {space.dimension}')
    import matplotlib.pyplot as plt
    from filmstar.plotting.plotter import update_plot
    try:
        figure = update_plot(projection, np.linspace(0, 1, 65))
        if args.save:
            figure.savefig(args.save, dpi=180)
            print(f'Saved {args.save}')
        else:
            plt.show(block=True)
    finally:
        plt.close('all')


if __name__ == '__main__':
    main()
