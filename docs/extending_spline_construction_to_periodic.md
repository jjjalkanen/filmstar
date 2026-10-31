**The same kernel construction works with periodic boundary conditions:** identify opposite boundary edges and add matching equations between their polynomial pieces, accounting for the translation used in the identification.

For example, take the rectangular fundamental domain

$\Omega =[0,{L}_{x}]\times [0,{L}_{y}],$

with periods ${L}_{x},{L}_{y}>0$. Identifying opposite sides gives a torus.

Suppose \\(p\_{\\mathrm L}(x,y)\\) and \\(p\_{\\mathrm R}(x,y)\\) are the polynomial pieces adjacent to corresponding portions of the left and right boundaries. For periodic ${C}^{r}$ matching, impose

\\\[ \\boxed{ p\_{\\mathrm R}(x+L\_x,y)-p\_{\\mathrm L}(x,y) \=x^{r+1}q\_x(x,y). } \\\]

Here ${q}_{x}$ is an unknown polynomial cofactor. Translating the right piece brings both pieces into coordinates where the identified boundary is $x=0$.

This equation enforces

\\\[ \\partial\_x^a\\partial\_y^b p\_{\\mathrm R}(L\_x,y) \= \\partial\_x^a\\partial\_y^b p\_{\\mathrm L}(0,y), \\qquad a+b\\le r. \\\]

Likewise, for corresponding bottom and top pieces,

\\\[ \\boxed{ p\_{\\mathrm T}(x,y+L\_y)-p\_{\\mathrm B}(x,y) \=y^{r+1}q\_y(x,y). } \\\]

The algorithm changes only in its mesh assembly:

1. Choose a fundamental domain and its polynomial cell decomposition.  
2. Pair opposite boundary segments under the periodic translations. Refine the decomposition if their subdivisions do not match.  
3. Add the usual divisibility equations for interior interfaces.  
4. Add the translated matching equations for the paired boundary segments.  
5. Expand coefficients and compute the resulting matrix kernel.

The translations are fixed, so these remain **linear equations in the unknown polynomial coefficients**. For pieces of total degree at most $d$, the displayed seam cofactors have degree at most $d-r-1$. Include all seam pairings, including those incident to periodically identified corners.

Homogenization also works: the translation becomes

\\\[ X\\longmapsto X+L\_xZ, \\\]

and the vertical seam equation becomes

\\\[ P\_{\\mathrm R}(X+L\_xZ,Y,Z)-P\_{\\mathrm L}(X,Y,Z) \=X^{r+1}Q\_x(X,Y,Z). \\\]

The qualification concerns the paper’s **dimension formulas**. Those use assumptions about a contractible planar domain and its boundary. A torus has different topology, so those formulas require a new homology calculation. The direct matrix-kernel computation already includes the periodic constraints and needs no such dimension formula.

