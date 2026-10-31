Here are the 2D coordinates and folding instructions for the **Tetrahedron** and **Cube** as single-piece fold-out nets.

### **1\. Tetrahedron (4 faces)**

We use an edge length of 2\. The height of an equilateral triangle is $h=\sqrt{3}\approx 1.732$. The net is a large triangle composed of 4 smaller triangles.  
**Corners (2D Coordinates):**

> * **C0:** (0, 0\)  
> * **C1:** (2, 0\)  
> * **C2:** (4, 0\)  
> * **C3:** (1, 1.732)  
> * **C4:** (3, 1.732)  
> * **C5:** (2, 3.464)

**Edges (Connects Corners):**

> * **e1:** C0-C1 \| **e2:** C1-C2 \| **e3:** C0-C3  
> * **e4:** C1-C3 \| **e5:** C1-C4 \| **e6:** C2-C4  
> * **e7:** C3-C4 \| **e8:** C3-C5 \| **e9:** C4-C5

**Faces (Bounded by Edges):**

> * **F1 (Center):** Bounded by e4, e5, e7 (Corners: C1, C4, C3)  
> * **F2 (Left):** Bounded by e1, e3, e4 (Corners: C0, C1, C3)  
> * **F3 (Right):** Bounded by e2, e6, e5 (Corners: C1, C2, C4)  
> * **F4 (Top):** Bounded by e7, e8, e9 (Corners: C3, C4, C5)

**Gluing Pairs (Outer Contour Edges to Join):** To fold into a 3D tetrahedron, join these outer edges:

> 1. **e3** (C0-C3) joins with **e8** (C3-C5)  
> 2. **e1** (C0-C1) joins with **e2** (C1-C2)  
> 3. **e6** (C2-C4) joins with **e9** (C4-C5)

### **2\. Cube (6 faces)**

We use a standard "cross" net with an edge length of 1\.  
**Corners (2D Coordinates):**

> * **C0:** (1, 0\) \| **C1:** (2, 0\)  
> * **C2:** (1, 1\) \| **C3:** (2, 1\)  
> * **C4:** (0, 1\) \| **C5:** (3, 1\)  
> * **C6:** (0, 2\) \| **C7:** (1, 2\)  
> * **C8:** (2, 2\) \| **C9:** (3, 2\)  
> * **C10:** (1, 3)\| **C11:** (2, 3\)  
> * **C12:** (1, 4)\| **C13:** (2, 4\)

**Edges (Connects Corners):**

> * *Verticals:* **v1** (C0-C2), **v2** (C1-C3), **v3** (C4-C6), **v4** (C2-C7), **v5** (C3-C8), **v6** (C5-C9), **v7** (C7-C10), **v8** (C8-C11), **v9** (C10-C12), **v10** (C11-C13)  
> * *Horizontals:* **h1** (C0-C1), **h2** (C4-C2), **h3** (C2-C3), **h4** (C3-C5), **h5** (C6-C7), **h6** (C7-C8), **h7** (C8-C9), **h8** (C10-C11), **h9** (C12-C13)

**Faces (Bounded by Edges):**

> * **F1 (Bottom):** h1, h3, v1, v2 (Corners: C0, C1, C3, C2)  
> * **F2 (Center-Front):** h3, h6, v4, v5 (Corners: C2, C3, C8, C7)  
> * **F3 (Left):** h2, h5, v3, v4 (Corners: C4, C2, C7, C6)  
> * **F4 (Right):** h4, h7, v5, v6 (Corners: C3, C5, C9, C8)  
> * **F5 (Center-Top):** h6, h8, v7, v8 (Corners: C7, C8, C11, C10)  
> * **F6 (Back):** h8, h9, v9, v10 (Corners: C10, C11, C13, C12)

**Gluing Pairs (Outer Contour Edges to Join):**

> 1. **v1** (C0-C2) $\leftrightarrow$ **h2** (C4-C2)  
> 2. **v2** (C1-C3) $\leftrightarrow$ **h4** (C3-C5)  
> 3. **h1** (C0-C1) $\leftrightarrow$ **h9** (C12-C13) *(Top of back face meets bottom of bottom face)*  
> 4. **v3** (C4-C6) $\leftrightarrow$ **v9** (C10-C12)  
> 5. **h5** (C6-C7) $\leftrightarrow$ **v7** (C7-C10)  
> 6. **v6** (C5-C9) $\leftrightarrow$ **v10** (C11-C13)  
> 7. **h7** (C8-C9) $\leftrightarrow$ **v8** (C8-C11)

### **3\. Octahedron (8 faces)**

Here is the complete single-piece blueprint for the **Octahedron** (8 faces).

To ensure the folded structure is solid and not just connected by fragile corner points, this layout uses a staggered chain of 8 equilateral triangles. We use an edge length of 2 and a triangle height of $h=\sqrt{3}\approx 1.732$.

**Corners (2D Coordinates)**

* **C0:** (3, \-1.732)  
* **C1:** (0, 0\)  
* **C2:** (2, 0\)  
* **C3:** (4, 0\)  
* **C4:** (6, 0\)  
* **C5:** (-1, 1.732)  
* **C6:** (1, 1.732)  
* **C7:** (3, 1.732)  
* **C8:** (5, 1.732)  
* **C9:** (2, 3.464)

**Edges (Connects Corners)**

* *Internal Creases (Fold lines):* **c1** (C1-C6), **c2** (C2-C6), **c3** (C2-C7), **c4** (C3-C7), **c5** (C3-C8), **c6** (C2-C3), **c7** (C6-C7)  
* *Outer Boundaries (Cut lines):* **e1** (C1-C5), **e2** (C5-C6), **e3** (C6-C9), **e4** (C7-C9), **e5** (C7-C8), **e6** (C4-C8), **e7** (C3-C4), **e8** (C0-C3), **e9** (C0-C2), **e10** (C1-C2)

**Faces (Bounded by Edges)**

* **F1:** Bounded by e10, c2, c1 (Corners: C1, C2, C6)  
* **F2:** Bounded by c3, c7, c2 (Corners: C2, C7, C6)  
* **F3:** Bounded by c6, c4, c3 (Corners: C2, C3, C7)  
* **F4:** Bounded by c5, e5, c4 (Corners: C3, C8, C7)  
* **F5:** Bounded by c7, e4, e3 (Corners: C6, C7, C9)  
* **F6:** Bounded by e9, e8, c6 (Corners: C0, C2, C3)  
* **F7:** Bounded by c1, e2, e1 (Corners: C1, C6, C5)  
* **F8:** Bounded by e7, e6, c5 (Corners: C3, C4, C8)

**Gluing Pairs (Outer Contour Edges to Join)** To fold into the 3D octahedron, join these outer edges together:

1. **e2** (C5-C6) joins with **e3** (C6-C9) *(Closes the top-left gap)*  
2. **e5** (C7-C8) joins with **e4** (C7-C9) *(Closes the top-right gap)*  
3. **e10** (C1-C2) joins with **e9** (C0-C2) *(Closes the bottom-left gap)*  
4. **e7** (C3-C4) joins with **e8** (C0-C3) *(Closes the bottom-right gap)*  
5. **e1** (C1-C5) joins with **e6** (C4-C8) *(Wraps the far-left edge to the far-right edge)*

**4\. Dodecahedron (12 faces)**

Here is the blueprint for a dodecahedron net with a strict edge length of 1\.

We map this using the golden ratio $\phi$ and the exact vertical increments for a unit pentagon.

### **Exact Constants**

* **Golden Ratio: $\phi =\frac{1+\sqrt{5}}{2}$**  
* **Upper height increment: ${s}_{1}=\sin (7{2}^{\circ })=\frac{\sqrt{2+\phi }}{2}$**  
* **Lower height increment: ${s}_{2}=\sin (3{6}^{\circ })=\frac{\sqrt{3-\phi }}{2}$**

### **1\. Flower A Corners (Bottom Bowl)**

Flower A is symmetric across the y-axis. Here are its 20 unique 2D coordinates:

* **Center Pentagon (F1):**  
  * **A0: $(0,{s}_{1}+{s}_{2})$**  
  * **A1: $(-\frac{\phi }{2},{s}_{1})$**  
  * **A2: $(-\frac{1}{2},0)$**  
  * **A3: $(\frac{1}{2},0)$**  
  * **A4: $(\frac{\phi }{2},{s}_{1})$**  
* **Left Petals (F2 & F3):**  
  * **A5: $(-\phi ,{s}_{1}+{s}_{2})$**  
  * **A6: $(-\frac{\phi +1}{2},2{s}_{1}+{s}_{2})$**  
  * **A7: $(-\frac{\phi -1}{2},2{s}_{1}+{s}_{2})$**  
  * **A8: $(-\frac{\phi +1}{2},-{s}_{2})$**  
  * **A9: $(-\frac{2\phi +1}{2},0)$**  
  * **A10: $(-\frac{\phi +2}{2},{s}_{1})$**  
* **Bottom & Right Petals (F4, F5, F6):**  
  * **A11: $(-\frac{\phi }{2},-{s}_{1})$**  
  * **A12: $(0,-{s}_{1}-{s}_{2})$**  
  * **A13: $(\frac{\phi }{2},-{s}_{1})$**  
  * **A14: $(\frac{\phi +1}{2},-{s}_{2})$**  
  * **A15: $(\frac{2\phi +1}{2},0)$**  
  * **A16: $(\frac{\phi +2}{2},{s}_{1})$**  
  * **A17: $(\phi ,{s}_{1}+{s}_{2})$**  
  * **A18: $(\frac{\phi +1}{2},2{s}_{1}+{s}_{2})$**  
  * **A19: $(\frac{\phi -1}{2},2{s}_{1}+{s}_{2})$**

### **2\. Flower B Corners (Top Bowl)**

Flower B is drawn as the mirror of Flower A across the line containing A19-A18 (at $y=2{s}_{1}+{s}_{2}$). The two flowers have separate vertex IDs, including B18 and B19, and are connected by the rim identifications below.

The formula for every point ${B}_{k}$ given ${A}_{k}=({x}_{k},{y}_{k})$ is:${B}_{k}=({x}_{k},4{s}_{1}+2{s}_{2}-{y}_{k})$

*B18 and B19 coincide with A18 and A19 in the flat drawing, but are distinct mesh vertices. Coinciding drawing coordinates do not create an uncut hinge.*

### **3\. Faces (12 Pentagons)**

All faces are perfectly regular with edge length 1\.

**Flower A:**

* **F1 (Center): A2, A3, A4, A0, A1**  
* **F2 (Top-Left): A0, A7, A6, A5, A1**  
* **F3 (Bottom-Left): A1, A10, A9, A8, A2**  
* **F4 (Bottom): A2, A11, A12, A13, A3**  
* **F5 (Bottom-Right): A3, A14, A15, A16, A4**  
* **F6 (Top-Right): A4, A17, A18, A19, A0**

**Flower B:**

* **F7 to F12: Mirror the exact same vertex sequences using the B-coordinates (e.g., F7 Center is B2, B3, B4, B0, B1).**

### **4\. Gluing Pairs**

To fold the 3D surface, glue the outer contour edges as follows:

**Stage 1: Close Flower A into a Bowl (Pinching the gaps between petals)**

1. **A0-A7 joins A0-A19**  
2. **A4-A17 joins A4-A16**  
3. **A3-A14 joins A3-A13**  
4. **A2-A11 joins A2-A8**  
5. **A1-A10 joins A1-A5**

**Stage 2: Close Flower B into a bowl**

6. B0-B7 joins B0-B19
7. B4-B17 joins B4-B16
8. B3-B14 joins B3-B13
9. B2-B11 joins B2-B8
10. B1-B10 joins B1-B5

**Stage 3: Join the two ten-edge rims with a one-edge offset**

11. A18-A19 joins B17-B18
12. A6-A7 joins B18-B19
13. A5-A6 joins B6-B7
14. A10-A9 joins B5-B6
15. A8-A9 joins B10-B9
16. A11-A12 joins B8-B9
17. A12-A13 joins B11-B12
18. A14-A15 joins B12-B13
19. A15-A16 joins B14-B15
20. A17-A18 joins B15-B16

The offset pairs each petal tip with a two-face rim vertex of the other bowl. After gluing, the mesh has 20 vertices, each incident to three edges and three pentagons.
