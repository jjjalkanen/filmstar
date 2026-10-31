# Polyhedron cutout splines

The four JSON files `examples/meshes/tetrahedron.json`, `examples/meshes/cube.json`,
`examples/meshes/octahedron.json`, and `examples/meshes/dodecahedron.json` encode the flat nets in
`periodic_test_cutouts.md`. Each face lists its ordered corners. Every cut
boundary edge has one partner in `edge_identifications`; the loader derives
the orientation-preserving rigid map by reversing the boundary direction of
the partner. Uncut shared edges use the existing interior-edge constraints.

For a glued edge, the solver imposes
`p_target(R X + a) - p_source(X) = b(X)^(smooth+1) q(X)`, where `R` is the
derived rotation and `b=0` defines the source edge. Thus C¹ matching is in
the source face's local edge directions. At an identified vertex, the solver
also takes two distinct boundary tangents of each incident face, crosses
their lifted 3D tangents, and compares the resulting normals after rotation
into one chart. The vertex equations are applied within the edge-constraint
nullspace. In the checked examples they do not reduce its dimension.

| Net | C⁰ degree 1 | C¹ degree 3 |
| --- | ---: | ---: |
| Tetrahedron | 4 | 4 |
| Cube | 4 | 4 |
| Octahedron | 6 | 5 |
| Dodecahedron | 4 | — |

The dodecahedron's C¹ degree 1 space has dimension 1, so it contains only
constant height functions. Higher degree C¹ spaces on that net have not been
measured. Exact nullspace calculations on glued nets stop after 55 seconds.
The corrected rim pairing has 20 quotient vertices of valence 3 and is also
the parameter mesh for the spatial dodecahedron example.

To produce a plot, run for example:

```sh
./run.sh --mesh examples/meshes/tetrahedron.json --deg 3 --smooth 1 \
  --save plots/tetrahedron_c1.png
```

Generated PNGs in `plots/` show the spline on the flat cutout coordinates. Glued
edges need not appear adjacent in that view; their values and derivatives
agree after applying the stored rigid map.
