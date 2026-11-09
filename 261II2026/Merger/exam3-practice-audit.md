# Exam 3 practice-material check

The question numbers were checked against [e33.pdf](e33.pdf), and every
answer in the original [e32.pdf](e32.pdf) was independently recomputed.
Both original PDFs are preserved. The packet now includes the editable
[corrected answer key](e32-corrected.tex)
([PDF](e32-corrected.pdf)).

## Concept-map references

| Topic | Practice questions in e33.pdf |
| --- | --- |
| Curl, divergence, conservativity | 1 |
| Potentials and FTLI | 2–5 |
| Green's theorem: circulation and planar flux | 6–13 |
| Scalar surface integrals, area, mass | 14–22 |
| Vector surface integrals / flux | 23–28 |
| Stokes on surfaces with boundary | 29–33 |
| Stokes on closed surfaces | 34–35 |
| Surface parametrizations | 36(a–h) |
| Operator types | 37 |
| Physical curl / divergence | 38 / 39 |
| Theorem statements and method selection | 40–41 |

The original diagram used a different numbering scheme. Its inferred,
cropped ranges have been replaced with these verified references.
The two editable map pages are [line integrals](exam3-map-lines.tex) and
[surface integrals](exam3-map-surfaces.tex).

## Answer-key findings

- **1(d) was incomplete.** It previously stated only “conservative.”
  The completed answer is
  \[
  \operatorname{curl}\vec F=(0,0,0),\qquad
  \operatorname{div}\vec F=-xz\cos y,
  \]
  with potential \(xz\cos y+C\).
- **No incorrect numerical answer or parametrization was found.** All
  other formulas, bounds, potentials, surface setups, and values in
  questions 1–41 check out under the stated course conventions.
- **41(3) and 41(5) need a method-choice qualification.** Their spatial
  curves lie in tilted planes. The original key's direct-method selections
  are retained, but the corrected key also identifies Green's theorem as
  valid after rewriting the integral on the planar projection. The claim
  “not in the xy-plane, so Green cannot apply” was too broad.
  In 41(3), the reduced form is
  \((y+2)\,dx+(x+y)\,dy\), with planar curl zero.
  In 41(5), it is
  \((y-x)\,dx+(1-2x-y)\,dy\), with planar curl \(-3\).

The verified values of 41(1–6), although not requested there, are
\(2\pi,0,0,2,-3/2,0\).

## Minor question-wording conventions

These are explained in the corrected key; `e33.pdf` itself is unchanged.

- **1(b):** The field excludes \(z=0\); it is nonconservative on each
  half-space \(z>0\) and \(z<0\).
- **17:** Read the surface scalar function as \(f(x,y,z)=y^5\).
- **29:** The intended hemisphere includes its equatorial rim, using
  \(z\ge0\), for the standard Stokes boundary calculation.
- **36(a):** The given bounds describe a rectangle, despite the word
  “square.” The original parametrization correctly follows the bounds.
- **37:** The answer uses the course convention that \(\nabla\) means
  the gradient of a scalar field, rather than a general vector Jacobian.

## Layout

[exam3mats.tex](exam3mats.tex) now produces **14 letter-sized pages**:
two separate concept-map pages, the formula sheet, overview, corrected
four-page key, and the original practice questions. Every page has the
same 612-by-792-point letter MediaBox. The enlarged A3 canvas is gone;
the imported formula sheet and practice PDFs retain their letter layout.
