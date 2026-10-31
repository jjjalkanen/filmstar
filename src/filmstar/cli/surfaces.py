"""Solve or load spatial spline patches and plot them in three dimensions."""
import argparse
import json
from pathlib import Path

import numpy as np

from filmstar.spatial.mesh import load_surface_mesh
from filmstar.spatial.solver import SurfaceSolution, SurfaceSpace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--deg', required=True, help='polynomial degree or auto (C1: 5, 6, 7)')
    parser.add_argument('--smooth', type=int, choices=(0, 1), required=True)
    parser.add_argument('--save', type=Path, metavar='PNG')
    parser.add_argument('--solution', type=Path, help='load saved coefficients')
    parser.add_argument('--write-solution', type=Path, help='save solved coefficients')
    args = parser.parse_args()
    mesh = load_surface_mesh(args.mesh)
    saved = json.loads(args.solution.read_text()) if args.solution else None
    if saved is not None:
        degrees = [int(saved['degree'])]
        if saved['mesh'] != mesh.name or saved['smoothness'] != args.smooth:
            raise ValueError('Saved solution does not match mesh or smoothness')
    elif args.deg == 'auto':
        degrees = [5, 6, 7] if args.smooth == 1 else [1]
    else:
        degrees = [int(args.deg)]
    if args.deg != 'auto' and int(args.deg) not in degrees:
        raise ValueError('Degree differs from saved solution')
    failure = None
    for degree in degrees:
        space = SurfaceSpace(mesh, degree, args.smooth)
        if saved is not None:
            coefficients = np.asarray(saved['coefficients'], dtype=float)
            if coefficients.shape != (3, space.columns):
                raise ValueError('Saved coefficient shape does not match degree')
            parameters = (coefficients-space.particular)@space.nullspace
            solution = SurfaceSolution(space, coefficients, parameters)
            break
        if args.smooth == 0:
            solution = space.c0_surface()
            break
        try:
            solution = space.solve_c1(max_nfev=600)
            break
        except ValueError as exc:
            failure = exc
    else:
        raise ValueError(f'No verified C1 surface through degree 7: {failure}')
    if np.max(np.abs(space.c0_matrix@solution.coefficients.T-space.c0_targets)) > 1e-8:
        raise ValueError('Loaded or computed solution fails C0 and vertex constraints')
    if args.smooth == 1:
        edge_error, normal_error, minimum = space.c1_errors(solution.coefficients)
        if edge_error > 1e-6 or normal_error > 1e-5 or minimum < 1e-8:
            raise ValueError('Loaded or computed solution fails C1 verification')
        if solution.independent_edge_error() > 1e-6:
            raise ValueError('Loaded or computed solution fails independent edge checks')
        metrics = solution.geometry_metrics()
        if (metrics['minimum_normal_ratio'] < 1e-3
                or (metrics['minimum_outward_cosine'] is not None
                    and metrics['minimum_outward_cosine'] < 0.05)
                or (metrics['maximum_radius_deviation'] is not None
                    and metrics['maximum_radius_deviation'] > 0.3)):
            raise ValueError(f'Loaded or computed surface folds or collapses: {metrics}')
    if args.write_solution:
        args.write_solution.write_text(json.dumps({
            'mesh': mesh.name, 'degree': degree, 'smoothness': args.smooth,
            'coefficients': solution.coefficients.tolist()}, indent=2)+'\n')
    print(f'{mesh.name}: {len(mesh.parameter_mesh.faces)} faces, '
          f'{len(mesh.vertex_ids)} spatial vertices, degree {degree}, C^{args.smooth}, '
          f'C0 free parameters {space.c0_dimension}', flush=True)

    import matplotlib
    matplotlib.use('Agg' if args.save else 'qtagg')
    import matplotlib.pyplot as plt
    from filmstar.plotting.plotter import update_plot
    from filmstar.plotting.spatial import SurfaceProjection
    try:
        figure = update_plot(SurfaceProjection(solution), np.linspace(0, 1, 65))
        if args.save:
            figure.savefig(args.save, dpi=180)
            print(f'Saved {args.save}', flush=True)
        else:
            plt.show(block=True)
    finally:
        plt.close('all')


if __name__ == '__main__':
    main()
