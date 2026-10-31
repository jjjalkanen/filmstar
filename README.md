# Filmstar

Filmstar constructs polynomial spline spaces on planar meshes, including
periodic and glued boundaries, and vector-valued spline surfaces on spatial
meshes. Equation solving and visualization live in separate modules.

## Install

Use Python 3.10 or newer. From the repository root:

```sh
python3 -m venv venv
venv/bin/python -m pip install -e .
```

The shell wrappers use `venv/bin/python` when available and otherwise use
`python3`. They add `src/` to the import path, so they also work with the
project dependencies installed without an editable install.

## Run examples

```sh
./run.sh --mesh examples/meshes/paper.json --deg 3 --smooth 0 \
  --save plots/paper_c0_d3.png
./run_surface.sh --mesh examples/surface_meshes/dodecahedron.json --deg 6 \
  --smooth 1 --solution examples/surface_meshes/dodecahedron_c1_solution.json \
  --save plots/dodecahedron_surface_c1.png
```

Omit `--save` to open an interactive Qt plot. Saved spatial solution JSON files
are included with the example meshes. Generated PNGs are ignored by Git.

## Test

```sh
venv/bin/python -m unittest discover -s tests -v
```

Tests include exact planar spline dimensions, glued and periodic examples,
spatial surface verification, and plotting command shutdown. Some exact solves
take several minutes.

Mesh format notes and the construction background are in `docs/` and
`examples/meshes/FORMAT.txt`.
