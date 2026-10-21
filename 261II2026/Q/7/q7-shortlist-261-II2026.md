# Quiz 7 — shortlist and teaching notes

Quiz date: Friday, October 23, 2026. Instructor notes; updated October 21.
Choices C, F, and revised H are confirmed as Problems 1, 2, and 4.
The product-rule question is confirmed as Problem 3. No point assignments.

Final files: [quiz source](q7-261-II2026.tex),
[quiz PDF](q7-261-II2026.pdf), [solution source](q7-sol-261-II2026.tex),
and [solution PDF](q7-sol-261-II2026.pdf).

These twelve candidates cover Green's theorem, conservative fields in
three dimensions, curl, and divergence. Each has proposed student wording,
a checked answer, and a source or adaptation note. The direct homework
alignment is with [HW7](../../HW/789/HW7_withSpace.tex): curl/divergence
(1–8), potentials and path independence (9–13), and Green's theorem for
circulation and outward flux (14–21).
[HW8](../../HW/789/HW8_withSpace.tex) and
[HW9](../../HW/789/HW9_withSpace.tex) were also reviewed; their surface
integrals, Stokes' theorem, and divergence theorem are deliberately outside
this shortlist's scope.

## At a glance

| Choice | Main idea | Suggested use |
| --- | --- | --- |
| A | Curl and divergence at a point; rule out conservativity | Core calculation |
| B | Choose a parameter to make divergence zero everywhere | Short calculation with two solutions |
| C | Choose parameters to make a 3D field conservative; find a potential | Selected: Problem 1 |
| D | A complicated spatial path becomes an endpoint calculation | Strong core quiz option |
| E | Green's theorem over a triangle | Standard Green calculation |
| F | Circulation versus outward flux for the same field | Selected: Problem 2 |
| G | Close a semicircle and subtract a nonzero segment integral | First open-path class example |
| H | A rounded-corner open path; close and subtract | Selected: Problem 4; revised October 21 |
| I | An open circular arc with a nonconstant Green integrand | Harder computational alternative |
| J | Four alternating semicircles | Visual classroom challenge |
| K | Repeated curl/divergence operations | Optional computational extension |
| L | Zero curl on a punctured domain is not enough | Optional conceptual extension |

The selected questions retain C and F exactly in mathematical content.
H now uses the requested three-piece path instead of the half-annulus.
The unused shortlist entries remain available for classroom examples.

## Selected Problem 3 — product rule and an unused alternative

### Selected: a product rule combined with divergence of a curl

**Proposed student wording.** Let
\[
f(x,y,z)=x^2+y^2+z^2,\qquad \vec G(x,y,z)=(z^2,x^2,y^2).
\]
Compute \(\nabla\cdot\bigl(f(\nabla\times\vec G)\bigr)\).

**Instructor check.** The product rule and the identity
\(\nabla\cdot(\nabla\times\vec G)=0\) give
\[
\begin{aligned}
\nabla\times\vec G&=(2y,2z,2x),\\
\nabla\cdot\bigl(f(\nabla\times\vec G)\bigr)
&=\nabla f\cdot(\nabla\times\vec G)
  +f\nabla\cdot(\nabla\times\vec G)\\
&=(2x,2y,2z)\cdot(2y,2z,2x)
=\boxed{4(xy+yz+xz)}.
\end{aligned}
\]
This tests exactly the operations and product rule covered in class, with
one short polynomial for the final-answer box. It is an application, not
a request to prove a general identity. A nonconstant scalar multiple of
a curl need not be divergence-free. Direct expansion is also a valid
method and gives the same result.

### Alternative: a very guided Helmholtz-type decomposition

**Proposed student wording.** Find real constants \(a,b\) such that
\[
\vec F=(2x-4y,\ 4x+2y,\ 2z)
=-\nabla f+\nabla\times\vec G,
\]
where
\[
f=a(x^2+y^2+z^2),\qquad
\vec G=(0,0,b(x^2+y^2)).
\]

**Instructor check.**
\[
-\nabla f+\nabla\times\vec G
=(-2ax+2by,\ -2ay-2bx,\ -2az),
\]
so \(\boxed{(a,b)=(-1,-2)}\). The prescribed forms make the two
coefficients unique; unrestricted scalar and vector potentials would not
be unique. This is a direct decomposition identity, not an invocation of
a global theorem or conditions at infinity.

The product-rule question is now Problem 3 of the quiz. The decomposition
is retained only as an unused alternative or a guided classroom example.

## A. A local computation that rules out a global potential

### Proposed student wording

Let
\[
\vec F(x,y,z)=(xz^2,\ x^2y-z^3,\ 2xy+y^2z).
\]
Find \((\nabla\cdot\vec F)(1,1,1)\) and
\((\nabla\times\vec F)(1,1,1)\). Is \(\vec F\) conservative on
\(\mathbb R^3\)?

### Instructor check

\[
\nabla\cdot\vec F=x^2+y^2+z^2,\qquad
\nabla\times\vec F=(2x+2yz+3z^2,\ 2xz-2y,\ 2xy).
\]
The answers are **\(3\)**, **\((7,0,2)\)**, and **no**. A nonzero curl at
even one point rules out a potential on all of \(\mathbb R^3\).

Source: [UCR surface-integral list](../../SourcesUCR/Lista-7-MA-1003-Int-superficie-II-2021.pdf#page=2),
page 2, exercise 1.5. Only the vector field is borrowed; the surface-integral
task is replaced with operations already in HW7.

## B. Make the field divergence-free

### Proposed student wording

For a real constant \(\lambda\), let
\[
\vec F_\lambda(x,y,z)
=((\lambda^2+x)^2,\ -4\lambda y,\ -2xz).
\]
Find all values of \(\lambda\) for which
\(\nabla\cdot\vec F_\lambda=0\) everywhere in \(\mathbb R^3\).

### Instructor check

\[
\nabla\cdot\vec F_\lambda
=2(\lambda^2+x)-4\lambda-2x
=2\lambda(\lambda-2).
\]
Answer: **\(\lambda\in\{0,2\}\)**.

Source: supplied “Chapter 7: Gradient, Curl, & Divergence — Practice
Problems” image, exercise 3. “Divergence zero everywhere” replaces the
flow-language interpretation so no additional terminology is required.
This is a useful shorter option, rather than the quiz's hardest calculation.

## C. Which parameters produce a conservative field?

### Proposed student wording

For real constants \(a,b\), let
\[
\vec F_{a,b}(x,y,z)
=(2xy+az^2,\ bx^2+z,\ 2xz+y).
\]
Find all pairs \((a,b)\) for which \(\vec F_{a,b}\) is conservative.
For those values, find a potential function.

### Instructor check

\[
\nabla\times\vec F_{a,b}
=(0,\ 2(a-1)z,\ 2(b-1)x).
\]
It vanishes everywhere exactly when **\((a,b)=(1,1)\)**. A potential is
\[
\boxed{f(x,y,z)=x^2y+xz^2+yz+C}.
\]
Differentiating this expression recovers all three components. Any
potential differing by an additive constant receives the same credit.

Source: original design, aligned with HW7's conservative-field and
potential questions. It joins a curl computation to potential recovery
without lengthy algebra.

## D. A complicated path with a simple endpoint calculation

### Proposed student wording

Let
\[
\vec F(x,y,z)=(3y^2z+ye^x,\ 6xyz+e^x,\ 3xy^2).
\]
Find a potential function and evaluate
\(\int_C\vec F\cdot d\vec r\), where
\[
\vec r(t)=
\left(\sin(t/2)+\cos(2t)-1,\ t\cos t,\ t^2\sin t\right),
\qquad 0\le t\le\pi,
\]
with \(C\) traversed in the direction of increasing \(t\).

### Instructor check

\[
f(x,y,z)=3xy^2z+ye^x+C,\qquad \nabla f=\vec F.
\]
The endpoints are \(\vec r(0)=(0,0,0)\) and
\(\vec r(\pi)=(1,-\pi,0)\). Therefore
\[
\boxed{\int_C\vec F\cdot d\vec r=f(1,-\pi,0)-f(0,0,0)=-\pi e}.
\]
All additive constants in the potential are acceptable. The field and
potential are smooth on all of \(\mathbb R^3\).

Source: [UCR line-integral list](../../SourcesUCR/Lista-6-MA-1003-Int-línea-II-2021.pdf#page=2),
page 2, exercise 1.8. The field and path are retained, and the potential
request is explicit. The first component contains \(\cos(2t)\), not
\(\cos^2t\). This closely parallels HW7's endpoint-method problems.

## E. Green's theorem on a triangle

### Proposed student wording

Use Green's theorem to evaluate
\[
\oint_C -y^2\,dx+x^2\,dy,
\]
where \(C\) is the counterclockwise boundary of the triangle with vertices
\((0,0)\), \((1,2)\), and \((0,2)\).

### Instructor check

The region has bounds
\[
\left\{
\begin{aligned}
&0\le y\le2,\\
&0\le x\le y/2.
\end{aligned}
\right.
\]
Since \(Q_x-P_y=2x+2y\),
\[
\oint_C -y^2\,dx+x^2\,dy
=\int_0^2\int_0^{y/2}(2x+2y)\,dx\,dy
=\boxed{\frac{10}{3}}.
\]

Source: triangle from [2025 Quiz 7, question 1](../../../261II2025/Q/7/q7-261-II2025.tex).
The old parametrization/bounds problem is rewritten as a Green's theorem
calculation. The variable inner bound remains important.

## F. Same field, same boundary: circulation and outward flux

### Proposed student wording

Let
\[
\vec F(x,y)=(y^3-2y,\ -4x+3y^2),
\]
and let \(C\) be the counterclockwise boundary of the rectangle
\(0\le x\le2\), \(0\le y\le1\). Find:

1. The circulation \(\oint_C\vec F\cdot\vec T\,ds\), where
   \(\vec T\) is the unit tangent in the direction of traversal.
2. The outward flux \(\oint_C\vec F\cdot\vec N\,ds\), where
   \(\vec N\) is the outward unit normal.

### Instructor check

\[
Q_x-P_y=-2-3y^2,\qquad P_x+Q_y=6y.
\]
Thus
\[
\text{circulation}
=\int_0^2\int_0^1(-2-3y^2)\,dy\,dx=\boxed{-6},
\qquad
\text{outward flux}
=\int_0^2\int_0^1 6y\,dy\,dx=\boxed{6}.
\]
The two requests would need separate labeled answer boxes in the final quiz.

Source: adapted from [2024 Quiz 8, question 2](../../../261II2024/Q/8/q8-261-II2024.tex),
with smaller powers, coefficients, and rectangle; orientation is explicit.
This matches the paired circulation/flux exercises in HW7.

## G. An open semicircle: the closing segment is not free

### Proposed student wording

Evaluate
\[
\int_C (y^2+x^3)\,dx+(2xy+x)\,dy,
\]
where \(C\) is the upper semicircle \(x^2+(y-1)^2=4\), \(y\ge1\),
traversed from \((2,1)\) to \((-2,1)\).

### Instructor check

Add the segment \(S:(-2,1)\to(2,1)\). The completed boundary is
counterclockwise, and \(Q_x-P_y=1\). The enclosed half-disk has area
\(2\pi\), while
\[
\int_S (y^2+x^3)\,dx+(2xy+x)\,dy
=\int_{-2}^2(1+x^3)\,dx=4.
\]
Consequently,
\[
\boxed{\int_C\vec F\cdot d\vec r=2\pi-4}.
\]

Source: original design inspired by the requested “add the missing segment”
method. This is my preferred first class example before H. The added
segment contributes a nonzero number, so students cannot omit it.

## H. A rounded-corner open path — selected Problem 4

### Proposed student wording

Use Green's theorem to evaluate \(\int_C\vec F\cdot d\vec r\), where
\[
\vec F(x,y)=(ye^{xy}-y,\ xe^{xy}+x+1).
\]
The path \(C\) follows, in order:

1. The straight segment from \((0,0)\) to \((3,0)\).
2. The counterclockwise quarter-circle \((x-2)^2+y^2=1\), from
   \((3,0)\) to \((2,1)\).
3. The straight segment from \((2,1)\) to \((0,1)\).

Reminder: \(\iint_R1\,dA=\operatorname{Area}(R)\).

### Instructor check

The specified center is \((2,0)\): this is the convex rounded corner,
not the other possible 90-degree arc centered at \((3,1)\).
Adding \(S:(0,1)\to(0,0)\) completes a counterclockwise boundary.
The enclosed region is a \(2\)-by-\(1\) rectangle plus a quarter of a
unit disk, so its area is \(2+\pi/4\). Equivalently, its bounds are
\[
\left\{
\begin{aligned}
&0\le y\le1,\\
&0\le x\le2+\sqrt{1-y^2}.
\end{aligned}
\right.
\]
Although the field looks unpleasant along the arc, its scalar curl is
\[
Q_x-P_y
=\bigl(e^{xy}+xy e^{xy}+1\bigr)
 -\bigl(e^{xy}+xy e^{xy}-1\bigr)=2.
\]
Thus
\[
\int_{C+S}\vec F\cdot d\vec r
=2\left(2+\frac\pi4\right)=4+\frac\pi2.
\]
On the added segment, \(x=0\), \(dx=0\), and \(Q(0,y)=1\), hence
\[
\int_S\vec F\cdot d\vec r=\int_1^0 1\,dy=-1.
\]
The requested integral is therefore
\[
\boxed{\int_C\vec F\cdot d\vec r
=\left(4+\frac\pi2\right)-(-1)=5+\frac\pi2}.
\]

**Why this field:** on the quarter-circle,
\(x=2+\cos t,\ y=\sin t\), the exponential becomes
\(e^{(2+\cos t)\sin t}\). Direct substitution is cumbersome, whereas
Green's theorem requires only cancellation, an area, and the constant
closing-segment integral. The field is smooth everywhere.

**Independent check:** the field can be written as
\(\nabla(e^{xy}+y)+(-y,x)\). The gradient part contributes
\((e^0+1)-(e^0+0)=1\) between the endpoints. The remaining part
contributes twice the area, \(4+\pi/2\), since its closing-segment
integral is zero. This gives the same total. Recognizing this splitting
is another mathematically valid shortcut; the student wording specifies
Green's theorem as the intended method.

Source: original field constructed for the instructor's October 21 path
variation, replacing the earlier half-annulus. The quiz includes a
black-and-white diagram of the actual open path, with direction arrows
and no added closing segment.

## I. A harder open arc with a nonconstant Green integrand

### Proposed student wording

Evaluate
\[
\int_C
\left(-\frac{y^3}{3}+2y-\sin(x+y)\right)\,dx
+
\left(\frac{x^3}{3}+2x-\sin(x+y)\right)\,dy,
\]
where \(C\) is the upper semicircle \(x^2+y^2=4x\), \(y\ge0\),
traversed from \((4,0)\) to \((0,0)\).

### Instructor check

Add \(S:(0,0)\to(4,0)\). The completed boundary is counterclockwise.
The trigonometric derivatives cancel, giving \(Q_x-P_y=x^2+y^2\).
Using \(x=2+r\cos\theta,\ y=r\sin\theta\), the half-disk bounds are
\[
\left\{
\begin{aligned}
&0\le\theta\le\pi,\\
&0\le r\le2.
\end{aligned}
\right.
\]
Hence
\[
\int_{C+S}\vec F\cdot d\vec r
=\int_0^\pi\int_0^2(r^2+4+4r\cos\theta)r\,dr\,d\theta
=12\pi.
\]
The closing segment contributes
\[
\int_S\vec F\cdot d\vec r
=\int_0^4-\sin x\,dx=\cos4-1.
\]
Answer: **\(12\pi+1-\cos4\)**, with angles in radians.

Source: [Moodle 009 question bank](../../SourcesUCR/Moodle/preguntas-II.S.2021.RRF.MA-1003.009-top-20260406-2207.html),
“Green 01,” question ID 11603673. The field is retained; the full circle
is replaced by its upper arc. Unlike H, this also requires a nonconstant
double integral, so it is the more computationally demanding choice.

## J. Four alternating semicircles — image-inspired classroom challenge

### Proposed student wording

All semicircles below lie above \(y=1\) and have center \((0,1)\).
Traverse the arcs in the listed order, joining the end of each arc to the
start of the next by a straight horizontal segment.

| Arc radius | Start | End |
| --- | --- | --- |
| \(4\) | \((4,1)\) | \((-4,1)\) |
| \(3\) | \((-3,1)\) | \((3,1)\) |
| \(2\) | \((2,1)\) | \((-2,1)\) |
| \(1\) | \((-1,1)\) | \((1,1)\) |

For the resulting open path \(C\), evaluate
\[
\int_C -2y\,dx+3x^2\,dy.
\]

### Instructor check

Add \(S:(1,1)\to(4,1)\). Its portion from \((2,1)\) to \((3,1)\)
cancels the original segment traversed from \((3,1)\) to \((2,1)\).
The remaining traversal is the sum of the positively oriented boundaries
of two upper half-annuli, with radii \(3\) to \(4\) and \(1\) to \(2\).
Apply Green's theorem separately to these two regions, rather than
treating the completed traversal as one simple closed curve.

Symmetry eliminates the \(6x\) term in \(Q_x-P_y=6x+2\), so
\[
\int_{C+S}\vec F\cdot d\vec r
=2\cdot\frac{\pi}{2}(16-9+4-1)=10\pi.
\]
Since \(\int_S\vec F\cdot d\vec r=\int_1^4-2\,dx=-6\),
the answer is **\(10\pi+6\)**.

Source: direct adaptation of the supplied “Path integrals vs. loop
integrals” image, translated upward by one unit. The original version
centered at the origin has answer \(10\pi\) and a zero closing-segment
contribution. It would make a useful precursor before this translated
version. Use a diagram with direction arrows if this becomes a handout;
this is a classroom challenge rather than my default short-quiz choice.

## K. Curl and divergence used more than once

### Proposed student wording

Compute:

1. \(\nabla\times(\nabla\times\vec F)\), where
   \(\vec F=(e^y,e^z,e^x)\).
2. \(\nabla\cdot[\nabla(\nabla\cdot\vec G)]\), where
   \(\vec G=(x^3y,y^3z,z^3x)\).

### Instructor check

First,
\[
\nabla\times\vec F=(-e^z,-e^x,-e^y),\qquad
\boxed{\nabla\times(\nabla\times\vec F)=(-e^y,-e^z,-e^x)}.
\]
Also, \(\nabla\cdot\vec G=3x^2y+3y^2z+3z^2x\); summing its three
second pure partial derivatives gives
\[
\boxed{\nabla\cdot[\nabla(\nabla\cdot\vec G)]=6(x+y+z)}.
\]

Source: [2025 Quiz 7, question 2](../../../261II2025/Q/7/q7-261-II2025.tex).
The second task also appears in the supplied Chapter 7 “More Problems”
image, exercise 7. This is an optional extension; either part could stand
alone as a shorter question. The first answer is a vector, the second
a scalar.

## L. Zero curl, but no global potential

### Proposed student wording

On \(D=\{(x,y,z):x^2+y^2>0\}\), let
\[
\vec F(x,y,z)=
\left(-\frac{y}{x^2+y^2},\ \frac{x}{x^2+y^2},\ 0\right).
\]
Find \(\nabla\times\vec F\) on \(D\) and
\(\int_C\vec F\cdot d\vec r\), where
\(\vec r(t)=(\cos t,\sin t,0)\), \(0\le t\le2\pi\), traversed in
increasing \(t\). Is \(\vec F\) conservative on \(D\)?

### Instructor check

The curl is **\((0,0,0)\)** everywhere on \(D\), but along \(C\),
\[
\vec F(\vec r(t))=(-\sin t,\cos t,0)=\vec r\,'(t).
\]
Thus the integral is **\(2\pi\)**, and the field is **not conservative**.
The domain excludes the entire \(z\)-axis. Green's theorem cannot be
applied to the filled unit disk because the field is undefined at its
center; a nonzero integral around a closed curve rules out a potential.

Source: three-dimensional extension of the supplied “Differentiability
matters” and “Orientation matters” examples. Use only if domain
hypotheses have been discussed; this is not a hidden prerequisite for
the core questions.

## Suggested teaching sequence and safeguards

Start with A or B, then C or D, then E or F. For the open-path technique,
use G as the first example and H as the harder follow-up. I is a different
computational version; J is the visual challenge closest to the supplied
image. K and L are optional extensions, not additional required topics.

For an open path \(C\) from \(A\) to \(B\), add a segment \(S\) from
\(B\) back to \(A\). If the completed boundary is positively oriented,
\[
\int_C P\,dx+Q\,dy
=\iint_D(Q_x-P_y)\,dA-\int_S P\,dx+Q\,dy.
\]
The orientation and the last subtraction both matter. G–J intentionally
have nonzero added-segment contributions; J additionally requires
boundary cancellation and decomposition into two regions.

Two cautions about the supplied/source material:

- The Chapter 6 practice image, exercise 3, has a sign typo in its flux
  formula: outward flux uses \(P_x+Q_y\), not \(P_x-Q_y\).
- Smoothness on the region matters, not just smoothness along the path.
  All core Green examples here use fields smooth on all of
  \(\mathbb R^2\). L explicitly illustrates the excluded-singularity
  issue; zero curl alone on an arbitrary domain is not enough.

C, F, the product-rule question, and revised H are now in the finalized
quiz. The solution key contains mathematical solutions only, without
grading notes. No homework or earlier quiz files have been changed.
