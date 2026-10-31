# Spatial polynomial spline surfaces

This is a separate extension of the planar height-function solver. A spatial
surface assigns a vector polynomial `(X(s,t), Y(s,t), Z(s,t))` to each face of
an existing flat parameter mesh. The 3D JSON supplies a `parameter_mesh` path
and exact `spatial_vertices` for every cutout vertex. Identified copies of a
vertex must have the same 3D position. Spatial edge curves are solved, not
fixed to the straight edges of a Platonic solid.

`SurfaceSpace(mesh, degree, 0)` constructs a C0 affine coefficient space:
`particular + nullspace @ parameters`. For C1, SymPy differentiates the face
monomials exactly. At each shared edge, the full inward transverse tangent
vectors are projected perpendicular to the common spatial edge tangent and
required to cancel. Their lengths are therefore matched in the flat face
charts. At each identified vertex, normals formed from two distinct incident
face edges are compared explicitly. These C1 equations are nonlinear in the
unknown coefficients; SciPy solves them using an analytic Jacobian. The C1
solution set is not presented as a linear nullspace.

For the circumscribed Platonic examples, the solver fixes nonzero edge
tangents at the vertices in their reference tangent planes. It uses a small
sphere-radius objective to select a
rounded surface when the supplied 3D vertices lie on a sphere. For other
closed meshes, it stays near the affine reference shape. It rejects surfaces
that fail independently sampled edge matching, collapse, or fold at sampled
points. This is a numerical check, not a proof against every self-intersection.

| Surface | Vertices | Faces | Saved C1 degree |
| --- | ---: | ---: | ---: |
| Tetrahedron | 4 | 4 | 5 |
| Cube | 8 | 6 | 5 |
| Octahedron | 6 | 8 | 5 |
| Dodecahedron | 20 | 12 | 6 |
| Torus | 9 | 9 | 6 |

The spatial dodecahedron uses `examples/meshes/dodecahedron.json` as its flat parameter
mesh. Its 20 quotient vertices all have valence 3.

The torus fixture uses a 3-by-3 periodic square parameter mesh. Its saved C1
surface uses cubic Hermite coordinate curves around each periodic direction;
their product has total degree 6. It verifies continuity without a spherical
reference. The general nonlinear search may still need a better initial
surface for other toroidal meshes, so use the saved solution for this example.

To plot a saved surface quickly:

```sh
./run_surface.sh --mesh examples/surface_meshes/dodecahedron.json --deg 6 --smooth 1 \
  --solution examples/surface_meshes/dodecahedron_c1_solution.json \
  --save plots/dodecahedron_surface_c1.png
```

Omit `--solution` to solve afresh. `--deg auto` tries degrees 5, 6, and 7 for
C1. Use `--smooth 0 --deg 1` to plot the faceted polyhedron. The shell wrapper
uses the project virtual environment, limits saved runs to 180 seconds, and
keeps numerical linear algebra single-threaded.
