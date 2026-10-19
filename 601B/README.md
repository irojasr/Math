# The Three Dictionaries of Toric Geometry

Notes for Renzo's Fall 2026 course at Colorado State University.

## Building the book

From the project directory, run:

```sh
pdflatex main.tex
makeindex main.idx
pdflatex main.tex
pdflatex main.tex
```

Alternatively, a standard `latexmk -pdf main.tex` workflow may be used.
The project requires a reasonably complete LaTeX installation with TikZ,
`tcolorbox`, `cleveref`, and the other packages listed in `preamble.tex`.

## Source layout

- `main.tex`: complete-book driver
- `preamble.tex`: packages, layout, theorem styles, and pedagogical boxes
- `macros.tex`: shared mathematical notation
- `frontmatter/`: title page, preface, and notation guide
- `chapters/`: Chapters 0 through 10
- `frontpageimage.png`: cover artwork

The master file divides Chapters 0--6 into Part I, *Basic Concepts*,
and Chapters 7--10 into Part II, *(More) Advanced Topics*. The
`makeindex` step generates the curated subject index printed at the end.

References are currently given in prose at the end of each chapter. The
project therefore does not require a `references.bib` file in this release.
