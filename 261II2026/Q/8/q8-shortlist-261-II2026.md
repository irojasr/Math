# Quiz 8 — shortlist

Quiz date: Friday, October 30, 2026. Instructor copy with checked answers.
Choices A--D are selected as Problems 1--4. No point assignments.

**The selected questions are A, B, C, and D:** one surface-area calculation, one
surface-mass calculation, one question on normal directions, and one on
the geometry of a parameter domain. Together they use all five families
of parametrizations discussed in class. The numerical answers have been
simplified to **\(42\pi\)** and **\(35\pi\)**, and the y-axis
parametrizations use the corrected outward-normal convention. E--H remain
unused alternatives.

Final files: [quiz source](q8-261-II2026.tex),
[quiz PDF](q8-261-II2026.pdf), [solution source](q8-sol-261-II2026.tex),
and [solution PDF](q8-sol-261-II2026.pdf).

| Choice | Topic | Main decision students must make |
| --- | --- | --- |
| **A** | **Area of a saddle-shaped sheet** | Find the surface-area factor and use a circular parameter region |
| **B** | **Mass of a conical lampshade** | Rotate a curve about the correct axis and substitute the density |
| **C** | **Which tangent vector comes first?** | Choose cross-product orders for outward or upward normals |
| **D** | **Two parametrizations of a triangular sail** | Distinguish a triangle from a parallelogram in parameter space |
| E | Area above a triangle | Include the changing width of the base region |
| F | Mass of a spherical cap | Translate the cutting plane into a bound on \(\phi\) |
| G | Mass of a cylinder with a sloping top | Use a height bound that depends on the angular parameter |
| H | Same rings, different connections | Recognize the surface produced by a ruled parametrization |

The calculations use scalar surface integrals and elementary integration.
The conceptual questions have short, identifiable answers suitable for
boxes; no proof or extended explanation is requested.

## A. Surface area of a saddle-shaped sheet — selected Problem 1

**Proposed student wording**

Find the area of the portion of the surface \(z=xy\) satisfying
\(x^2+y^2\le15\).

**Instructor check**

For \(\vec r(x,y)=(x,y,xy)\),
\[
\vec r_x\times\vec r_y=(-y,-x,1),\qquad
\|\vec r_x\times\vec r_y\|=\sqrt{1+x^2+y^2}.
\]
The disk is described in polar coordinates by
\[
\left\{
\begin{aligned}
&0\le\theta\le2\pi,\\
&0\le r\le\sqrt{15}.
\end{aligned}
\right.
\]
Thus
\[
A=\int_0^{2\pi}\int_0^{\sqrt{15}}r\sqrt{1+r^2}\,dr\,d\theta
=\frac{2\pi}{3}\bigl(16^{3/2}-1\bigr)
=\boxed{42\pi}.
\]
Equivalently, students may start directly with
\(\vec r(r,\theta)=(r\cos\theta,r\sin\theta,r^2\cos\theta\sin\theta)\).

**Why I like it.** The surface is curved, but its area factor becomes
radial. The radius \(\sqrt{15}\) makes the final answer an integer multiple
of \(\pi\). This builds on
the plane-over-a-disk example without repeating its constant area factor.

**Source.** Adapted from the
[UCR surface-integral list](../../SourcesUCR/Lista-7-MA-1003-Int-superficie-II-2021.pdf#page=2),
page 2, exercise 1.1; the original cylinder has radius \(\sqrt8\).

## B. A lampshade that gets heavier toward the top — selected Problem 2

**Proposed student wording**

A thin lampshade is formed by rotating the line segment
\(y=\frac34x\), \(2\le x\le4\), about the \(y\)-axis. Its surface density is
\(\sigma(x,y,z)=y\). Find the mass of the curved surface, without end caps.

**Instructor check**

Using the corrected convention, one parametrization is
\[
\vec r(t,\theta)=(t\cos\theta,\tfrac34t,-t\sin\theta),\qquad
\left\{
\begin{aligned}
&0\le\theta\le2\pi,\\
&2\le t\le4.
\end{aligned}
\right.
\]
Then
\[
\vec r_\theta\times\vec r_t
=(\tfrac34t\cos\theta,-t,-\tfrac34t\sin\theta),\qquad
\|\vec r_\theta\times\vec r_t\|=\tfrac54t,
\qquad \sigma(\vec r)=\tfrac34t.
\]
Therefore
\[
M=\int_0^{2\pi}\int_2^4\frac{15}{16}t^2\,dt\,d\theta
=\boxed{35\pi}.
\]
The cross product displayed here is outward: its dot product with the
radial vector \((t\cos\theta,0,-t\sin\theta)\) is \(3t^2/4>0\).
Reversing its order leaves the scalar area element and mass unchanged.
Other parametrizations covering the same lateral surface once give the
same mass. The slope \(3/4\) gives a \(3\)-\(4\)-\(5\) triangle and
eliminates the square root from the area factor.

**Why I like it.** Students must distinguish the planar generating
coordinate from the spatial coordinates after rotation. The nonconstant,
positive density makes this a mass problem with a short polynomial integral.

**Source.** Original variation aligned with
[HW8](../../HW/789/HW8_withSpace.tex), exercises 5, 7, 9, and 10.

## C. Which tangent vector comes first? — selected Problem 3

**Proposed student wording**

For each parametrization, write the order of the tangent-vector cross
product that gives the requested normal direction. Use the right-hand
rule; no cross-product calculation is required.

1. **Outward** from the ellipsoid:
   \[
   \vec r(\theta,\phi)
   =(2\sin\phi\cos\theta,3\sin\phi\sin\theta,\cos\phi),
   \qquad 0\le\theta<2\pi,\quad0<\phi<\pi.
   \]
2. **Outward** from the \(y\)-axis:
   \[
   \vec q(t,\theta)=(\cos\theta,t,-\sin\theta),
   \qquad0\le\theta<2\pi,\quad0<t<2.
   \]
3. **Upward** on the graph:
   \[
   \vec p(u,v)=(v,u,u^2+v^2).
   \]

Write each answer in the form \(\vec r_\theta\times\vec r_\phi\), using
the appropriate parametrization and order.

**Instructor check**

The three answers are
\[
\boxed{\vec r_\phi\times\vec r_\theta},\qquad
\boxed{\vec q_\theta\times\vec q_t},\qquad
\boxed{\vec p_v\times\vec p_u}.
\]
For verification,
\[
\begin{aligned}
\vec r_\phi\times\vec r_\theta
&=(3\sin^2\phi\cos\theta,2\sin^2\phi\sin\theta,
6\sin\phi\cos\phi),\\
\vec q_\theta\times\vec q_t&=(\cos\theta,0,-\sin\theta),\\
\vec p_v\times\vec p_u&=(-2v,-2u,1).
\end{aligned}
\]
The first is a positive multiple of the outward normal to
\(x^2/4+y^2/9+z^2=1\); the second points away from the cylinder's axis;
the third has positive \(z\)-component. The poles are excluded because
the ellipsoid's angular parametrization degenerates there.

**Why I like it.** The cylinder is around the \(y\)-axis, and the graph's
parameters are deliberately exchanged. Students must read the actual
parametrization rather than automatically cross its first two labels.
The quiz uses three small, separately labeled answer boxes, following
the format of the previous orientation question.

**Source.** Adapted from
[Quiz 9 (2025), Question 3](../../../261II2025/Q/9/q9-261-II2025.tex).
That is the normal-orientation question you remembered. There is also a
sphere/tangent-direction exercise in
[Quiz 9 (2024)](../../../261II2024/Q/9/q9-261-II2024.tex).

## D. Two students, one triangular sail — selected Problem 4

**Proposed student wording**

A triangular sail has vertices
\[
P=(0,0,1),\qquad Q=(2,0,1),\qquad R=(0,3,4).
\]
Two proposed parametrizations are
\[
\begin{aligned}
\vec a(s,t)&=P+s(Q-P)+t(R-P),
&&0\le s\le1,\quad0\le t\le1,\\
\vec b(u,v)&=(1-v)P+v\bigl((1-u)Q+uR\bigr),
&&0\le u\le1,\quad0\le v\le1.
\end{aligned}
\]

(a) Which covers exactly the triangular sail: **\(\vec a\) only,
\(\vec b\) only, both, or neither**?

(b) What additional inequality on \(s,t\) makes \(\vec a\) cover exactly
the triangular sail?

**Instructor check**

Answers: **\(\vec b\) only**, and **\(s+t\le1\)**.

The square domain for \(\vec a\) produces the parallelogram with fourth
vertex \(Q+R-P=(2,3,4)\). For \(\vec b\), the point
\((1-u)Q+uR\) runs along the edge \(QR\), and varying \(v\) draws the
segment from \(P\) to that point. Its image is exactly the triangular sail.
Equivalently, the corrected bounds for \(\vec a\) are
\[
\left\{
\begin{aligned}
&0\le s\le1,\\
&0\le t\le1-s.
\end{aligned}
\right.
\]
The ruled parametrization collapses the edge \(v=0\) to \(P\); this does
not change its image or prevent its use for surface integration.

**Why I like it.** It joins the plane and ruled-surface constructions and
tests whether students understand what the bounds do. It also follows
naturally from your explanation of triangular domains. Two short boxes
suffice: a choice and an inequality.

**Source.** Original design using the two constructions taught in class;
related to the triangular base in [HW8](../../HW/789/HW8_withSpace.tex),
exercise 1(d).

## E. A curved sheet over a triangle — alternative area question

**Proposed student wording**

Find the area of the portion of \(z=x^2/2\) above the triangular region
\(0\le y\le x\le1\).

**Instructor check**

Here \(dS=\sqrt{1+x^2}\,dx\,dy\). The bounds are
\[
\left\{
\begin{aligned}
&0\le x\le1,\\
&0\le y\le x.
\end{aligned}
\right.
\]
Consequently
\[
A=\int_0^1\int_0^x\sqrt{1+x^2}\,dy\,dx
=\int_0^1x\sqrt{1+x^2}\,dx
=\boxed{\frac{2\sqrt2-1}{3}}.
\]
A unit-square parametrization is
\(\vec r(u,v)=(v,uv,v^2/2)\), with \(0\le u,v\le1\);
its area factor is \(v\sqrt{1+v^2}\).

**Why choose it instead of A?** It tests a triangular base rather than
polar coordinates, and the changing slice width supplies the factor
needed for a simple substitution.

**Source.** Original variation on the graph-over-a-triangle construction
in HW8, exercise 1(d), with a different surface and base triangle.

## F. Mass of a spherical cap — alternative mass question

**Proposed student wording**

Find the mass of the portion of the sphere \(x^2+y^2+z^2=4\) with
\(z\ge1\), if its surface density is \(\sigma(x,y,z)=z\).
Include only the curved surface.

**Instructor check**

Use
\[
\vec r(\theta,\phi)
=(2\sin\phi\cos\theta,2\sin\phi\sin\theta,2\cos\phi),
\qquad
\left\{
\begin{aligned}
&0\le\theta\le2\pi,\\
&0\le\phi\le\pi/3.
\end{aligned}
\right.
\]
Then \(dS=4\sin\phi\,d\phi\,d\theta\), giving
\[
M=\int_0^{2\pi}\int_0^{\pi/3}
8\cos\phi\sin\phi\,d\phi\,d\theta
=\boxed{6\pi}.
\]

**Why choose it instead of B?** The main challenge is finding the polar
angle bound. This gives more weight to spherical parametrizations while
keeping the mass calculation short.

**Source.** Adapted from HW8, exercise 14, with a changed cutting plane
and a nonconstant density.

## G. A cylindrical wall with a sloping top — harder mass alternative

**Proposed student wording**

Find the mass of the portion of the cylinder \(x^2+z^2=4\) satisfying
\(0\le y\le3+x\), with surface density \(\sigma(x,y,z)=y\).
Include only the curved wall.

**Instructor check**

Take
\[
\vec r(\theta,t)=(2\cos\theta,t,2\sin\theta),\qquad
\left\{
\begin{aligned}
&0\le\theta\le2\pi,\\
&0\le t\le3+2\cos\theta.
\end{aligned}
\right.
\]
The height bound is always positive. Since the area factor is \(2\),
\[
M=\int_0^{2\pi}\int_0^{3+2\cos\theta}2t\,dt\,d\theta
=\int_0^{2\pi}(3+2\cos\theta)^2\,d\theta
=\boxed{22\pi}.
\]

**Why choose it?** It tests a nonrectangular parameter domain and a
cylinder around the \(y\)-axis. It is longer than B
or F. If selected, the fact
\(\int_0^{2\pi}\cos^2\theta\,d\theta=\pi\) could be supplied so that
trigonometric integration does not dominate the question.

**Source.** Adapted from HW8, exercise 8, changing the cylinder axis,
cutting plane, and task from area setup to mass evaluation.

## H. Same rings, different connections — creative conceptual alternative

**Proposed student wording**

Two unit rings lie in the planes \(z=0\) and \(z=2\). Straight strings
join the points
\[
\vec r_1(s)=(\cos s,\sin s,0),\qquad
\vec r_2(s)=(-\cos s,-\sin s,2)
\]
with matching values of \(s\). The strings form the ruled surface
\[
\vec R(s,t)=(1-t)\vec r_1(s)+t\vec r_2(s),\qquad
0\le s\le2\pi,\quad0\le t\le1.
\]

(a) Which shape do they form: **a cylinder, a conical frustum, or two
cones meeting at a vertex**?

(b) Find the point at which all the strings meet.

**Instructor check**

\[
\vec R(s,t)=((1-2t)\cos s,(1-2t)\sin s,2t).
\]
At \(t=1/2\), every string passes through **\((0,0,1)\)**. Eliminating
the parameters gives
\[
x^2+y^2=(1-z)^2,\qquad0\le z\le2.
\]
The shape is **two cones meeting at a vertex**. Joining matching angular
points without the minus signs would instead produce a cylinder.
The vertex is singular, but this question only asks students to identify
the image, not construct a normal there.

**Why I like it.** The two boundary rings do not determine the ruled
surface by themselves: the correspondence between points matters.
This is a visual application of the ruling formula, with almost no
calculation. It could replace D if you want a less familiar conceptual
question. A final-quiz sketch should show the two rings and the pairing
of one pair of endpoints, leaving the complete shape for students to infer.

**Source.** Original variation inspired by the ruled-surface construction
in [Quiz 9 (2025), Question 2](../../../261II2025/Q/9/q9-261-II2025.tex).

## Instructor note: the corrected y-axis revolution convention

The selected convention is
\[
\vec f(t,\theta)=(x(t)\cos\theta,y(t),-x(t)\sin\theta).
\]
Its tangent vectors are \(\vec f_t=M(\theta)\vec a'(t)\) and
\(\vec f_\theta=M'(\theta)\vec a(t)\), where
\[
M(\theta)=
\begin{pmatrix}
\cos\theta&0\\
0&1\\
-\sin\theta&0
\end{pmatrix},\qquad
\vec a(t)=\begin{pmatrix}x(t)\\y(t)\end{pmatrix}.
\]
In particular,
\[
\vec f_\theta\times\vec f_t
=(xy'\cos\theta,-xx',-xy'\sin\theta).
\]
Its dot product with the radial vector
\((x\cos\theta,0,-x\sin\theta)\) is \(x^2y'\). Thus, when
\(x>0\) and \(y'>0\), this order points **outward**. The direction in
which the generating curve is traversed matters. For scalar area and
mass, either cross-product order has the same norm.

## Sources reviewed and scope

- [HW8](../../HW/789/HW8_withSpace.tex), especially 1(b, d, f), 5–10,
  and 14. The shortlist builds on the reported class examples without
  repeating them verbatim.
- [Quiz 8 (2025)](../../../261II2025/Q/8/q8-261-II2025.tex): graph
  area setup and a spherical peel parametrization.
- [Quiz 9 (2025)](../../../261II2025/Q/9/q9-261-II2025.tex): ruled
  surfaces and choosing cross-product orders.
- [Quiz 9 (2024)](../../../261II2024/Q/9/q9-261-II2024.tex): sphere
  area and tangent-vector orientation.
- [UCR surface-integral list](../../SourcesUCR/Lista-7-MA-1003-Int-superficie-II-2021.pdf):
  exercise 1.1 supplies the saddle-area idea. The later flux and theorem
  questions are outside this quiz's present scope.

All final values, cross-product directions, and parameter bounds above
have been checked. The selected questions and their worked solutions are
in the linked final quiz and solution files above.
