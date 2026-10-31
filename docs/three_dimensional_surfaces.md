We **can use the same general approach to construct polynomial surface patches in \\(\\mathbb R^3\\), with free coefficients remaining after the matching equations are solved.** The extra ingredient is specifying how the two-dimensional coordinates on adjacent patches relate across their shared edge.

The planar construction already contains a simple instance of this: a scalar spline $f(x,y)$ produces the surface

\\\[ \\mathbf X(x,y)=(x,y,f(x,y)). \\\]

For a general surface, make all three coordinates unknown.

For each face ${\sigma }_{i}$ of an abstract polyhedral mesh, choose a two-dimensional parameter domain ${D}_{i}$, with coordinates $u,v$, and write

\\\[ \\boxed{ \\mathbf X\_i(u,v) \= \\sum\_{a+b\\le d}\\mathbf c\_{i,a,b}u^av^b, \\qquad \\mathbf c\_{i,a,b}\\in\\mathbb R^3. } \\\]

Here \\(\\mathbf X\_i\\) maps the parameter face into space, $d$ is the degree bound, and the vectors \\(\\mathbf c\_{i,a,b}\\) are the unknowns. The mesh specifies which patches meet; its initial polyhedral shape need not be the resulting surface.

Suppose two patches share an edge, represented in both parameter domains by $(t,0)$, with corresponding points using the same $t$. If its prescribed spatial curve is \\(\\boldsymbol\\gamma\_e(t)\\), impose

\\\[ \\boxed{ \\mathbf X\_i(t,0)=\\mathbf X\_j(t,0) \=\\boldsymbol\\gamma\_e(t). } \\\]

These are three polynomial identities per equality. They fix the boundary curves while leaving interior shape parameters available, provided the degree and other constraints allow them.

For example, on a triangular parameter domain, let its three boundary lines have equations \\(\\ell\_1,\\ell\_2,\\ell\_3\\). Once a patch \\(\\mathbf X\_i^0\\) with the desired boundary is found, every modification

\\\[ \\mathbf X\_i \= \\mathbf X\_i^0+\\ell\_1\\ell\_2\\ell\_3\\,\\mathbf Q\_i \\\]

preserves that boundary, where \\(\\mathbf Q\_i\\) is any vector polynomial within the degree budget. Thus the boundary curves generally do not determine the patch interior.

**For continuity alone, that is essentially the construction. For derivative matching, we also need gluing data.**

In the planar case, adjacent pieces already share the coordinates $x,y$. On separately parameterized surface patches, their transverse derivatives need to be compared through a specified coordinate transition. A first-order matching condition can take the form

\\\[ \\boxed{ \\partial\_v\\mathbf X\_i(t,0) \= a\_e(t)\\,\\partial\_v\\mathbf X\_j(t,0) \+ b\_e(t)\\,\\partial\_u\\mathbf X\_j(t,0). } \\\]

The functions ${a}_{e},{b}_{e}$ describe how the transverse parameter direction on one patch corresponds to directions on the other. With compatible, fixed polynomial gluing data, this is linear in the surface coefficients. Rational gluing data can also be used, clearing denominators that are nonzero along the edge.

This is the established framework of **geometrically continuous splines**: polynomial patches connected through transition maps. Algebraic versions explicitly extend the Billera–Schenck–Stillman construction and provide basis algorithms. [arxiv.org](https://arxiv.org/abs/2305.09096?utm_source=chatgpt.com)

There are three concrete qualifications:

* **Spatial edge curves alone do not specify derivative matching.** They specify the common boundary and its tangent, but leave the transverse directions undetermined. Fixing the gluing data preserves a linear system; solving for those data together with the patches generally introduces nonlinear equations.  
* **Prescribed edges must permit the requested smoothness at vertices.** At a regular ${C}^{1}$ surface point, every incident edge tangent must lie in one tangent plane. Three linearly independent prescribed edge tangents cannot meet at such a point. For example, preserving all three perpendicular cube edges at a corner is incompatible with making that corner a regular smooth surface point.  
* **Matching equations produce a surface map, not automatically an embedding.** Regularity requires\\\[ \\partial\_u\\mathbf X\_i\\times\\partial\_v\\mathbf X\_i\\ne0, \\\]and an embedding additionally excludes self-intersections.

So, with parameter faces, boundary curves, degree, and compatible gluing data fixed, the procedure remains: **assemble coefficient equations, solve an affine linear system, and use its nullspace to parameterize the remaining surfaces.** Merely replacing planar edge equations by equations of space curves omits the patch coordinates and derivative identifications needed to define that system.

**The gluing data specify how parameter directions and derivative scales correspond across an edge.** For a smooth join, both patches must have the same physical tangent plane there. An actual angle between their tangent planes would produce a crease; that is also possible, but it is a different matching requirement.

**At a corner, compatibility involves every incident patch**, not necessarily three. The coordinate transitions must agree after going around the vertex. At first order, if ${T}_{1},\ldots ,{T}_{k}$ are the $2\times 2$ derivative matrices of successive transitions, the consistency condition is

${T}_{k}\cdots {T}_{2}{T}_{1}=I.$

Higher smoothness requires the corresponding higher-order transition compatibility. Prescribed spatial edge tangents must also fit the desired tangent plane at the vertex.

With compatible gluing data fixed, yes: choose polynomial degrees, write the boundary and derivative matching identities, and solve for the patch coefficients. The result has the form

\\\[ \\boxed{ \\mathbf c=\\mathbf c\_0+\\sum\_j\\lambda\_j\\mathbf n\_j, } \\\]

where \\(\\mathbf c\_0\\) is one solution, the \\(\\mathbf n\_j\\) form a nullspace basis, and the ${\lambda }_{j}$ are the remaining free shape parameters.

The qualification is existence: a chosen degree and prescribed boundaries may yield no solution, or only geometrically degenerate solutions. The equation construction still works; obtaining a regular embedded surface additionally requires selecting parameters that avoid collapsed patches and self-intersections.

