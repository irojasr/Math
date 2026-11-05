# Quiz 9 — shortlist

Quiz date: Friday, November 6, 2026. Instructor copy with checked answers
and source notes. **A, E, C, and G are selected as Problems 1--4**, in
that order. No point values are assigned.

The selected quiz has two conceptual questions, a direct flux calculation,
and a Stokes calculation. B, D, F, and H remain unused alternatives for
classroom examples.

Final files: [quiz source](q9-261-II2026.tex),
[quiz PDF](q9-261-II2026.pdf), [solution source](q9-sol-261-II2026.tex),
and [solution PDF](q9-sol-261-II2026.pdf).

| Choice | Topic | Expected answer |
| --- | --- | --- |
| **A — Problem 1** | **Which side of Stokes do we compute?** | Line integral of the field; surface integral of its curl |
| B | An enclosed region and a surface with no boundary | Braced inequalities; empty boundary |
| **C — Problem 3** | **Flux through a tilted triangular sheet** | \(3\) |
| D | Stokes on a half-ellipsoid with a complicated field | \(-6\pi\) |
| **E — Problem 2** | **Two boundary circles and their directions** | Outer CCW, inner CW |
| F | Direct flux through a cylindrical wall | \(30\pi\) |
| **G — Problem 4** | **Stokes turns large powers into a constant curl** | \(6\) |
| H | Stokes on a tilted circular boundary | \(24\pi\) |

Throughout, \(\vec N\) denotes the specified **unit normal**, \(\vec T\)
the unit tangent in the direction of traversal, and \(d\sigma\) surface
area. All calculation answers are integers or integer multiples of
\(\pi\). These questions use vector surface integrals and Stokes;
the divergence theorem is not required.

## A. Which side of Stokes do we compute? — selected Problem 1

### Proposed student wording

In each part, use Stokes' theorem to replace the requested integral by
the **other side of the theorem**. State the **type of integral** and
the **vector field being integrated**. Do not compute derivatives or
evaluate either integral.

**(a)** Let \(\vec F=(y,z,x)\), and let \(S\) be the outward-oriented
upper hemisphere
\[
x^2+y^2+z^2=4,\qquad z\ge0.
\]
The requested integral is
\[
\iint_S (\operatorname{curl}\vec F)\cdot\vec N\,d\sigma.
\]
Does its Stokes replacement integrate \(\vec F\) or
\(\operatorname{curl}\vec F\), along a curve or over a surface?

**(b)** Let \(\vec G=(z,x,y)\), and let \(C\) traverse the triangle
\[
(2,0,0)\longrightarrow(0,2,0)\longrightarrow(0,0,2)
\longrightarrow(2,0,0).
\]
The requested integral is
\[
\oint_C\vec G\cdot\vec T\,ds.
\]
Does its Stokes replacement integrate \(\vec G\) or
\(\operatorname{curl}\vec G\), along a curve or over a surface?

### Expected answers and check

- **(a) A line integral of \(\vec F\).** The equality is
  \[
  \iint_S (\operatorname{curl}\vec F)\cdot\vec N\,d\sigma
  =\oint_{\partial S}\vec F\cdot\vec T\,ds.
  \]
  Here \(\partial S\) is \(x^2+y^2=4,\ z=0\), traversed
  counterclockwise as viewed from the positive \(z\)-axis.
- **(b) A surface integral of \(\operatorname{curl}\vec G\).** Taking
  \(S_T\) to be the flat triangular sheet with the given vertices,
  \[
  \oint_C\vec G\cdot\vec T\,ds
  =\iint_{S_T}(\operatorname{curl}\vec G)\cdot\vec N\,d\sigma.
  \]
  The plane is \(x+y+z=2\), and the matching normal is upward.

The student answers can be just the two bold phrases. The boundary and
surface descriptions above verify the replacements; they are not extra
requests in the proposed question. For an optional instructor check,
the values would be \(-4\pi\) and \(6\), respectively.

**Why this wording matters.** Directly evaluating the original integral
is mathematically valid. Asking for the **other side of Stokes** tests
the switch you want, with an unambiguous answer. In (b), integrating
\(\vec G\) over a surface would be the wrong replacement: the integrand
must involve its curl.

**Sources.** Adapted from [HW9](../../HW/789/HW9_withSpace.tex), problem
13 (hemisphere radius changed), and [Baohua's prequiz](BaohuaQ9Review.tex),
problem 5 (triangle scaled by 2). The hemisphere includes its rim so
Stokes applies to the stated surface with that boundary.

## B. An enclosed region and a surface with no boundary — unused conceptual alternative

### Proposed student wording

**(a)** Write bounds for the planar region \(D\) enclosed by
\(x=y^2-1\) and \(x=3\), with \(y\) as the outer variable.

**(b)** Let \(S\) be the **entire spherical surface**
\(x^2+y^2+z^2=4\). Give its boundary curves \(\partial S\), in the
sense used in Stokes' theorem. Write \(\varnothing\) if there are none.

### Expected answers and check

\[
\boxed{
\left\{
\begin{aligned}
&-2\le y\le2,\\
&y^2-1\le x\le3.
\end{aligned}
\right.
}
\qquad
\boxed{\partial S=\varnothing}.
\]
The two planar curves meet when \(y^2-1=3\), giving \(y=\pm2\).
For each such \(y\), the region runs from the parabola to the vertical
line. The full sphere has no edge curves, unlike a hemisphere.

The sphere is itself the boundary of a **solid ball**, but has no boundary
curves **as a surface**. This distinction makes the intended answer clear.
If the same sphere appears in a curl-flux integral of a smooth field,
Stokes therefore gives zero; that conclusion does not apply to arbitrary
flux integrals of \(\vec F\).

**Source and prior precedent.** The planar region is a shifted version
of your proposed parabola-and-line example. You have used the same
boundary/region question style before:
[Quiz 7 (2025), Question 1](../../../261II2025/Q/7/q7-261-II2025.tex)
asks for a triangle's boundary parametrization and inequalities for its
enclosed region. The sphere case aligns with HW9, problems 12 and 15,
and Baohua's conceptual problem 2.

For the final quiz, use a box large enough for the two-line bounds and a
small separate box for the empty boundary. A small planar sketch could
accompany the solution. "Enclosed region" is deliberate: its non-strict
bounds include the boundary; the strict topological interior would use
strict inequalities.

## C. Direct flux through a tilted triangular sheet — selected Problem 3

### Proposed student wording

Find the flux of \(\vec F=(y,x,y+z)\) through the portion of
\(2x+y+z=2\) in the first octant, oriented **upward**:
\[
\iint_S\vec F\cdot\vec N\,d\sigma.
\]

### Expected answer and check

Use \(\vec r(x,y)=(x,y,2-2x-y)\). The upward oriented area vector is
\[
\vec r_x\times\vec r_y=(2,1,1).
\]
The base triangle has bounds
\[
\left\{
\begin{aligned}
&0\le x\le1,\\
&0\le y\le2-2x.
\end{aligned}
\right.
\]
On the surface,
\[
\vec F(\vec r)=(y,x,2-2x),\qquad
\vec F(\vec r)\cdot(2,1,1)=2+2y-x.
\]
Therefore
\[
\begin{aligned}
\iint_S\vec F\cdot\vec N\,d\sigma
&=\int_0^1\int_0^{2-2x}(2+2y-x)\,dy\,dx\\
&=\int_0^1(8-14x+6x^2)\,dx
=\boxed{3}.
\end{aligned}
\]
The cross product already includes the area factor. No additional factor
\(\sqrt6\) is needed after taking this dot product.

**Why I like it.** It joins a familiar plane parametrization, the correct
normal direction, a triangular domain, and a genuine flux dot product.
The final integral is a short polynomial calculation.

**Source.** Adapted from the
[UCR surface-integral list](../../SourcesUCR/Lista-7-MA-1003-Int-superficie-II-2021.pdf#page=3),
page 3, exercise 1.13. The field and surface are retained; the original
asks for flux of the **curl** via Stokes. This version asks directly for
flux of \(\vec F\), so its answer is different.

## D. A complicated field, a simple boundary — unused Stokes alternative

### Proposed student wording

Let \(S\) be the outward-oriented half-ellipsoid
\[
9x^2+4y^2+36z^2=36,\qquad z\ge0,
\]
and let
\[
\vec F(x,y,z)=
\bigl(y,\ x^2,\ z(x^2-y^3)^7\cos(e^{xyz})\bigr).
\]
Use Stokes' theorem to find
\[
\iint_S(\operatorname{curl}\vec F)\cdot\vec N\,d\sigma.
\]
You may use
\(\int_0^{2\pi}\sin^2t\,dt=\pi\) and
\(\int_0^{2\pi}\cos^3t\,dt=0\).

### Expected answer and check

The boundary is the ellipse
\[
9x^2+4y^2=36,\qquad z=0,
\]
counterclockwise when viewed from above. A matching parametrization is
\[
\vec r(t)=(2\cos t,3\sin t,0),\qquad0\le t\le2\pi.
\]
Stokes replaces the surface integral by the line integral of \(\vec F\):
\[
\begin{aligned}
\vec F(\vec r(t))&=(3\sin t,4\cos^2t,0),\\
\vec r'(t)&=(-2\sin t,3\cos t,0),\\
\vec F(\vec r(t))\cdot\vec r'(t)
&=-6\sin^2t+12\cos^3t.
\end{aligned}
\]
Consequently,
\[
\iint_S(\operatorname{curl}\vec F)\cdot\vec N\,d\sigma
=\int_0^{2\pi}(-6\sin^2t+12\cos^3t)\,dt
=\boxed{-6\pi}.
\]
The complicated third component contributes nothing on this boundary;
no derivatives of it are needed. An outward normal on the upper half
induces the counterclockwise boundary direction used above.

**Why I like it.** It makes the benefit of the conceptual switch in A
immediately visible, while revisiting last week's ellipsoid parametrization
geometry. Students only need the ellipse at \(z=0\), not a parametrization
of the entire surface.

**Sources.** Adapted from
[Quiz 10 (2024), Question 1](../../../261II2024/Q/10/q10-261-II2024.tex)
and [Quiz 9 (2025), Question 4](../../../261II2025/Q/9/q9-261-II2025.tex).
The ellipsoid's horizontal semiaxes are interchanged from the 2024
version, and the third component follows the 2025 construction. Orientation
is stated explicitly.

## E. A boundary with two components — selected Problem 2

### Proposed student wording

Let \(S\) be the annular surface
\[
1\le x^2+y^2\le4,\qquad z=0,
\]
with upward normal. Give the equations of its boundary curves and the
induced direction of each, viewed from above.

### Expected answers and check

- **Outer boundary:** \(x^2+y^2=4,\ z=0\), **counterclockwise**.
- **Inner boundary:** \(x^2+y^2=1,\ z=0\), **clockwise**.

With the upward normal, the annular surface stays on the left when each
boundary component is traversed. It has two boundary circles, not just
the outer one.

**Source.** [Baohua's prequiz](BaohuaQ9Review.tex), problem 6(a).
Selected in place of B to emphasize induced orientation. The quiz uses
separate labeled boxes for the outer and inner boundaries.

## F. Direct flux through a cylindrical wall — alternative calculation

### Proposed student wording

Find the outward flux of \(\vec F=(xz,y,x^2)\) through the curved wall
\(x^2+y^2=4\), \(0\le z\le3\). Do not include end caps.

### Expected answer and check

With
\[
\vec r(\theta,z)=(2\cos\theta,2\sin\theta,z),\qquad
\left\{
\begin{aligned}
&0\le\theta\le2\pi,\\
&0\le z\le3,
\end{aligned}
\right.
\]
the outward oriented area vector is
\(\vec r_\theta\times\vec r_z=(2\cos\theta,2\sin\theta,0)\).
Hence
\[
\begin{aligned}
\iint_S\vec F\cdot\vec N\,d\sigma
&=\int_0^{2\pi}\int_0^3
(4z\cos^2\theta+4\sin^2\theta)\,dz\,d\theta\\
&=18\pi+12\pi=\boxed{30\pi}.
\end{aligned}
\]
Only the curved wall is included. The third field component does not
contribute because the wall's normal has zero \(z\)-component.

**Source.** [Baohua's prequiz](BaohuaQ9Review.tex), problem 4, with the
task changed from setup-only to evaluation. This is the strongest
alternative to C if you want a curved rather than planar flux surface.

## G. Large powers disappear under curl — selected Problem 4

### Proposed student wording

Let
\[
\vec F=(z-x^{2026},\ x-y^{2026},\ y-z^{2026}),
\]
and let \(C\) traverse
\[
(2,0,0)\longrightarrow(0,2,0)\longrightarrow(0,0,2)
\longrightarrow(2,0,0).
\]
Use Stokes' theorem to evaluate
\(\displaystyle\oint_C\vec F\cdot\vec T\,ds\).

### Expected answer and check

\[
\operatorname{curl}\vec F=(1,1,1).
\]
Choose the triangular sheet in the plane \(z=2-x-y\). The stated vertex
order corresponds to the upward oriented area vector \((1,1,1)\).
The projected bounds are
\[
\left\{
\begin{aligned}
&0\le x\le2,\\
&0\le y\le2-x.
\end{aligned}
\right.
\]
Thus
\[
\oint_C\vec F\cdot\vec T\,ds
=\int_0^2\int_0^{2-x}3\,dy\,dx
=\boxed{6}.
\]
The large powers disappear because only mixed-component derivatives
enter the curl. This is a concrete line-to-surface counterpart to D.

**Source.** Adapted from the
[Moodle .009 export](../../SourcesUCR/Moodle/preguntas-II.S.2021.RRF.MA-1003.009-top-20260406-2207.html),
"Circulación 1," question 13306988. The source has exponent 2021 and
parameterized intercepts \(2a\); here the intercepts are fixed at 2,
the exponent is 2026, and orientation is specified by vertex order.

## H. A tilted circular boundary — alternative Stokes calculation

### Proposed student wording

Let \(\vec F=(-y^3,x^3,-z^3)\), and let \(C\) be the intersection of
\[
x^2+y^2=4,\qquad z=2-x-y,
\]
traversed counterclockwise as viewed from above. Use Stokes' theorem to
evaluate \(\displaystyle\oint_C\vec F\cdot\vec T\,ds\).

### Expected answer and check

Take the upward-oriented planar sheet
\(\vec r(x,y)=(x,y,2-x-y)\) over the disk \(x^2+y^2\le4\).
Then
\[
\operatorname{curl}\vec F=(0,0,3(x^2+y^2)),\qquad
\vec r_x\times\vec r_y=(1,1,1).
\]
Therefore
\[
\oint_C\vec F\cdot\vec T\,ds
=\int_0^{2\pi}\int_0^2 3r^3\,dr\,d\theta
=\boxed{24\pi}.
\]
The spatial boundary is an ellipse in a tilted plane, with a circular
projection onto the \(xy\)-plane. The projected circle determines the
parameter region; it is not the boundary curve itself in three dimensions.

**Source.** Adapted from the
[UCR surface-integral list](../../SourcesUCR/Lista-7-MA-1003-Int-superficie-II-2021.pdf#page=3),
page 3, exercise 1.15. The cylinder radius and plane constant are changed,
and the orientation is specified as viewed from above. This is a little
more involved than G because the curl integrand is nonconstant.

## Reference notes

- [BaohuaQ9Review.tex](BaohuaQ9Review.tex) is an unchanged copy of the
  supplied `Prequiz9_Surface integraals of vector function+Stokes' Theorem.tex`.
  Its conceptual questions 1--3, setup questions 4--6, and computations
  8--9 are directly relevant to the current scope.
- [HW9](../../HW/789/HW9_withSpace.tex), problems 1--9, supplies vector
  flux examples; problems 10--18 supply Stokes examples, including closed
  surfaces, caps, and surfaces with more than one boundary component.
- Earlier geometry practice includes the triangle in
  [Quiz 7 (2025), Question 1](../../../261II2025/Q/7/q7-261-II2025.tex)
  and the parabolic/linear bounds in
  [Quiz 4 (2026)](../4/q4-261-II2026.tex).
- The UCR list and both Moodle exports were reviewed. Moodle's
  "Flujo sobre Superficies" (question 14214789 in the .008 export) is
  also a useful source for comparing surfaces and boundary components,
  but its embedded diagrams depend on remote image files; the candidates
  here give their geometry explicitly.
- [HW10](../../HW/10/HW10_withSpace.tex) and the remaining questions in
  the old Quiz 10 files move into the divergence theorem, so they were
  not used as requirements for these questions.

The selected questions and their worked solutions are in the final files
linked above. Unused entries remain available for classroom examples.
