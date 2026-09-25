# Quiz 4 — draft choices and alternatives

Instructor selection notes, revised September 23, 2026. The student draft is
[q4-261-II2026.tex](q4-261-II2026.tex), with a provisional date of September 25.
It retains four questions, choose exactly two, two pages, scratch space, and
unlabeled answer boxes at the ends of the problems. The matching solution key is
[q4-sol-261-II2026.tex](q4-sol-261-II2026.tex).

September 23 revision: Questions 1 and 3 have been swapped: Question 1 combines
two double integrals into one, and Question 3 splits one into two. Question 4
uses candidate B's direct triple-integral calculation. The quiz
does not ask students to change the order of a triple integral; those problems
are deferred until that topic has been covered.

## Questions currently in the draft

### 1. Combine two double integrals into one by reversing the order

Let \(f\) be continuous. Rewrite
\[
\int_0^1\int_0^{\sqrt{x}}f(x,y)\,dy\,dx
+\int_1^2\int_0^{2-x}f(x,y)\,dy\,dx
\]
as one iterated integral in the order \(dx\,dy\). Do not evaluate.

Source: original replacement, selected September 23 after the previous
\(\sin(y^2)\) order-reversal problem was used in class.

Expected answer:
\[
\int_0^1\int_{y^2}^{2-y}f(x,y)\,dx\,dy.
\]

The original pieces meet at \(x=1\). For each \(0\le y\le1\), their
horizontal intervals are \([y^2,1]\) and \([1,2-y]\), which join to form
\([y^2,2-y]\). This complements Question 3: there a reversal splits one
integral into two, while here it combines two integrals into one.

### 2. Mass setup over a triangular base with a curved top

A solid lies in the first octant, below \(z=4-y^2\), with \(x+y\le2\).
Its density is \(\rho(x,y,z)=1+z\). Set up a triple integral for its mass in
the order \(dz\,dy\,dx\). Do not evaluate.

Source: Moodle export
[section .008](../../SourcesUCR/Moodle/preguntas-II.S.2021.RRF.MA-1003.008-top-20260406-2208.html),
**VIT-V2**, question **11712500**, source comment at line 3185.
Adaptation: the original asks for volume; this version adds density and asks
only for the mass integral.

Expected answer:
\[
\int_0^2\int_0^{2-x}\int_0^{4-y^2}(1+z)\,dz\,dy\,dx.
\]

This is the requested bounds check. In Cartesian coordinates, the constant
limits \(0\le x,y\le2\), \(0\le z\le4\) describe a larger box. The triangular
base requires \(y\le2-x\), and the top requires \(z\le4-y^2\).

### 3. Double-order reversal with a split — selected candidate C

Let \(h\) be continuous. Rewrite
\[
\int_0^1\int_x^{2-x}h(x,y)\,dy\,dx
\]
in the order \(dx\,dy\). Do not evaluate.

Source: original variant of the linear-bound reversal practiced in class;
candidate **C** below, now selected to replace the exponential-triangle problem.

Expected answer:
\[
\int_0^1\int_0^y h(x,y)\,dx\,dy
+\int_1^2\int_0^{2-y}h(x,y)\,dx\,dy.
\]

The triangle has vertices \((0,0),(0,2),(1,1)\). Its horizontal slices change
at \(y=1\). An equivalent single integral uses \(0\le y\le2\) and
\(0\le x\le\min\{y,2-y\}=1-|y-1|\).

### 4. A short triple calculation with a small twist — selected candidate B

Let \(g\) be continuously differentiable, with \(g(0)=1\) and \(g(1)=4\).
Evaluate
\[
\int_0^1\int_0^1\int_0^{\sqrt{2x}}2g'(y)z\,dz\,dx\,dy.
\]

Source: Moodle [section .008](../../SourcesUCR/Moodle/preguntas-II.S.2021.RRF.MA-1003.008-top-20260406-2208.html),
**Jesus-SG-Integrales-triples-calculo-directo-P01-V01 (copiar)**,
question **11522301**, source comment at line 2343.
Adaptation: endpoint values are supplied to make the answer numerical.

Expected answer: \(3\). Integrating in \(z\) gives \(2xg'(y)\), and
integrating in \(x\) leaves
\(\int_0^1g'(y)\,dy=g(1)-g(0)=4-1=3\).
No change in the order of a triple integral is required.

## Other candidates

### A. Gentler triple-order reversal — deferred until covered in class

Let \(g\) be continuous. Rewrite
\[
\int_0^5\int_0^2\int_0^{\sqrt{4-y^2}}g(x,y,z)\,dx\,dy\,dz
\]
in the order \(dz\,dy\,dx\). Do not evaluate.

Source: Moodle [section .009](../../SourcesUCR/Moodle/preguntas-II.S.2021.RRF.MA-1003.009-top-20260406-2207.html),
**IOIT-V9**, question **11603684**, source comment at line 980.

Expected answer:
\[
\int_0^2\int_0^{\sqrt{4-x^2}}\int_0^5 g(x,y,z)\,dz\,dy\,dx.
\]

The cross section is a quarter-disk of radius 2 at every height. This candidate
is deferred until changes in the order of triple integrals have been covered;
polar coordinates are not needed.

### B. A short triple calculation with a small twist — now selected as Question 4

Let \(g\) be continuously differentiable, with \(g(0)=1\) and \(g(1)=4\).
Evaluate
\[
\int_0^1\int_0^1\int_0^{\sqrt{2x}}2g'(y)z\,dz\,dx\,dy.
\]

Source: Moodle [section .008](../../SourcesUCR/Moodle/preguntas-II.S.2021.RRF.MA-1003.008-top-20260406-2208.html),
**Jesus-SG-Integrales-triples-calculo-directo-P01-V01 (copiar)**,
question **11522301**, source comment at line 2343.
Adaptation: endpoint values are supplied to make the answer numerical.

Expected answer: \(3\). The first two integrations leave
\(\int_0^1g'(y)\,dy=g(1)-g(0)\).

Selected as Question 4 on September 23. It replaces the triple-order reversal
with a direct calculation and an application of the Fundamental Theorem of
Calculus.

### C. Double-order reversal with a split — now selected as Question 3

Let \(h\) be continuous. Rewrite
\[
\int_0^1\int_x^{2-x}h(x,y)\,dy\,dx
\]
in the order \(dx\,dy\). Do not evaluate.

Source: original variant of the linear-bound reversal practiced in class.

Expected answer:
\[
\int_0^1\int_0^y h(x,y)\,dx\,dy
+\int_1^2\int_0^{2-y}h(x,y)\,dx\,dy.
\]

The horizontal slices change at \(y=1\). This is now Question 3 in the quiz
and solution key, replacing the triangle computation from the 2025 quiz.

### D. A triangle from the 2024 quizzes

Let \(T\) be the triangle with vertices \((-2,2)\), \((-2,-2)\), and \((2,2)\).
Evaluate
\[
\iint_T 3(y-x)\,dA.
\]

Source: [2024 Quiz 5 solutions](../../../261II2024/Q/5/q5-sol-261-II2024.tex),
Problem **2**, beginning at line **221**. The original asks for a sketch,
the edges, and area integrals in both orders. This adaptation retains the
triangle and replaces those tasks with one weighted integral.

Expected answer: \(32\), using
\[
\int_{-2}^2\int_x^2 3(y-x)\,dy\,dx.
\]

This remains an alternative; candidate C is the selected Question 3.

## Notes on the supplied materials

- Baohua's A4 answer calls \(dz\,dy\,dx\) “\(z\) outermost.” In that notation,
  \(z\) is innermost. The draft names the differential order explicitly.
- The classroom example with inner limits \(\sqrt y\) to \(y\) is a signed
  integral: for \(0<y<1\), \(y<\sqrt y\). The usual positively oriented
  region integral has inner limits \(y\) to \(\sqrt y\).
- The draft uses Cartesian coordinates and the topics explicitly described
  in class. Polar problems from the review and Moodle remain available for
  a later quiz.
