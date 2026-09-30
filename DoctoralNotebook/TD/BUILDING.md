# Building the thesis

Run these commands from `DoctoralNotebook/TD/`. They use the existing
`make`, `latexmk`, pdfLaTeX, and BibTeX tools; no new LaTeX package is needed.

| Task | Command | PDF |
| --- | --- | --- |
| Complete thesis | `make thesis` (or `make`) | `build/td.pdf` |
| Introduction | `make introduction` | `build/introduction-preview.pdf` |
| Background | `make background` | `build/background-preview.pdf` |
| Main construction chapter | `make chapter1` | `build/chapter1-preview.pdf` |
| Remove all generated output, including PDFs | `make clean` | — |

For a complete rebuild before submission, run `make clean` followed by
`make thesis`. Check `build/td.log`, `build/td.blg`, and the resulting PDF.
The document's existing metadata and content TODOs still need to be completed.

## How previews work

A chapter contains `\chapter`, labels, and content, but no document class,
preamble, or document environment. Compile through `td.tex`; compiling a
chapter file directly will not supply those requirements.

The preview targets pass a chapter selector to `td.tex` using latexmk's
`-usepretex` option. The main applies `\includeonly` and uses exactly the same
class, packages, macros, fonts, margins, and theorem definitions as a complete
build. Only the selected body chapter is typeset. Preliminary pages are
omitted; the bibliography is included. The class's copyright-page counter is
initialized in preview mode so saved chapter counters can be restored.

The first preview after cleaning automatically runs a complete thesis build
to establish labels and chapter/page/counter checkpoints. Subsequent previews
reuse that snapshot without typesetting the other chapters. Each preview
copies the full-build auxiliary files into its own output directory and asks
latexmk to recompile the selected chapter and settle references and citations.
Preview output cannot overwrite the full thesis's saved checkpoints.

References to omitted chapters and their citation records come from the last
complete build. They resolve when that snapshot contains the referenced labels
and citations. Newly added labels or changes in omitted chapters need another
`make thesis`; otherwise references can be stale or undefined. Rebuild fully
after changing chapter order or counter/theorem definitions. A preview retains
the selected chapter's full-thesis number and starting page from the snapshot.
If its length changes, the bibliography's page numbering may reflect old
checkpoints from later omitted chapters until a complete rebuild.

## Bibliography order

The backend is BibTeX with `cite`, and the style is now `plain`: numeric
citations with entries sorted by author names (surname first), then year,
then title. The previous `unsrt` style used citation order; `\nocite{*}` made
the database order visible. Reordering `.bib` records is not the sorting
configuration for `plain`.

The supplied template describes `unsrt` as tested and permits other styles
(`docs/template-example.tex.txt`, bibliography discussion). The class's numeric
labels and bibliography layout are unchanged. All 30 current records, keys,
and fields in `bibiTD.bib` are preserved. A comparison of generated entries
under `plain` and `unsrt` found identical entry text, with only order changed.
Citation numbers change to match alphabetical order. The temporary
`\nocite{*}` still prints the whole database; remove it when the chapters have
actual citations. Cached citations from omitted chapters normally keep preview
reference numbers consistent with the full build's citation set.

## Output and editor use

All generated files live under ignored `build/`:

```text
build/
  td.pdf, td.aux, td.bbl, td.blg, td.log, td.toc, ...
  chapters/*.aux
  introduction-preview.pdf, background-preview.pdf, chapter1-preview.pdf
  previews/<chapter>/<chapter>-preview.*
  previews/<chapter>/chapters/*.aux
```

The per-preview auxiliary directories prevent LaTeX's identically named
chapter `.aux` files from colliding. Only output is separated; the chapter
sources remain together in `chapters/`. Latexmk also copies preview SyncTeX
data alongside the convenient PDFs in `build/`.

The local `.latexmkrc` also makes `latexmk td.tex` put a full build in `build/`.
Set an editor's root document to `td.tex` and use latexmk from this directory;
do not select `chapters/*.tex` as independent documents. Raw `pdflatex` does not
read `.latexmkrc`, so use the documented targets to keep output out of sources.
Relative `figures/` paths continue to resolve from the thesis directory.

## Why this approach

| Option | Assessment |
| --- | --- |
| Main document with `\includeonly` | Selected: native LaTeX, full formatting fidelity, cached labels/counters, clean chapter sources. |
| `subfiles` or `standalone` | Adds a package and chapter wrappers; no benefit for these existing `\include` chapters. |
| Separate preview driver | Would duplicate or extract the preamble and require more care to preserve chapter counters. |
| Editor root directive or latexmk command alone | Keeps the shared main, but a root directive alone does not select a chapter or isolate preview output. |
| Small build configuration | Selected alongside `\includeonly`: `.latexmkrc` controls output and a Makefile supplies readable commands and auxiliary-file snapshots. |

Verification covered the full thesis and all three previews, alphabetical
bibliography output, unique keys, and clean source directories. A separate
temporary fixture also tested real citations without `\nocite{*}`, references
to omitted chapters, theorem numbering, shared macros, and a PDF loaded from
`figures/`. Test content and assets were not added to the thesis chapters.
