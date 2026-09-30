Thesis project setup
===================

Build from DoctoralNotebook/TD/:

    make thesis        # build/td.pdf
    make introduction  # build/introduction-preview.pdf
    make background    # build/background-preview.pdf
    make chapter1      # build/chapter1-preview.pdf
    make clean         # remove generated build/ output, including PDFs

See BUILDING.md for the full workflow, output locations, bibliography sorting,
and preview limitations. The engine is pdfLaTeX; the bibliography processor is
BibTeX with plain (alphabetical, numeric citations). No Biber or biblatex is used.
The directory is self-contained; a normal TeX distribution and make are required.
The local .latexmkrc also sends `latexmk td.tex` output into build/. Chapter
previews use the same main preamble with includeonly and snapshots from the
last full build. Source chapters contain no generated auxiliary files.

Generated build/td.pdf is the full preview, not a source file. Local .gitignore
rules exclude build/ and legacy auxiliary outputs without ignoring figure PDFs.

Project files
-------------

Modified:
  td.tex                       Main document, metadata, matter boundaries,
                               chapter imports, bibliography, graphics paths.
  bibiTD.bib                   Consolidated and minimally repaired database.

Created:
  VarillyThesisMacros.tex       72 math commands and seven theorem environments.
  chapters/introduction.tex    Introduction; label chap:introduction.
  chapters/background.tex      Background; label chap:background.
  chapters/chapter1.tex        Pseudostable Maps and the Main Construction;
                               label chap:pseudostable-maps; title tentative.
  figures/.gitkeep             Keeps the empty figure directory in version control.
  docs/template-example.tex.txt Exact original td.tex, including all instructions,
                               sample prose, comments, and public-domain license.
                               Documentation only; not another build entry point.
  README.txt                   Build instructions and setup/audit report.
  .gitignore                   Generated-output exclusions local to this folder.
  .latexmkrc                   Shared pdfLaTeX defaults and build/ output.
  Makefile                     Full/preview/clean targets with isolated caches.
  BUILDING.md                  Current compilation and preview instructions.

thesis.cls is unchanged. No source outside TD was modified. The original
headers and DoctoralNotebook/bibiDoctoralNotebook.bib were not changed.

Architecture decision
---------------------

A: Keep td.tex plus one reduced mathematical header (selected).
   Directly recognizes the supplied university template; one build target;
   low conflict risk; portable folder; content and mathematical notation are
   separated while university formatting remains in the class. Reviewers can
   follow the metadata and front/main/back matter in one small main document.

B: obligatoryHeader.tex plus VarillyReduced.tex.
   Could be compatible, but moves the already-small template preamble behind
   another import. It adds maintenance and review indirection without a
   demonstrated portability or formatting benefit here.

C: mainTD.tex plus archived td.tex.
   Could preserve an example, but creates two apparent main documents and
   ambiguity for editors, build tools, and reviewers. Archiving the example
   as documentation instead preserves the instructions without that ambiguity.

Each chapter uses \include: chapters naturally start on new pages. The Makefile
now selects individual chapters through \includeonly in the main document.
After changing chapter structure, run a full build to refresh preview snapshots.
The main has no thesis prose. Chapter files contain only titles, unique labels,
and TODO comments, with no invented mathematical statements.

Original template and class
---------------------------

The supplied class identifies itself as the Colorado State University Thesis
class, dated 2020/05/07, based on book.cls. The source uses [doctor].
It prints DISSERTATION and Doctor of Philosophy for that option.

The original main was both an example and the package manual: five chapters
(Introduction, The Thesis Document Class, Figures and Tables, Creating
References and Citations, Formatting Tips and Tricks), bibliography, and a
License appendix. Its sample front matter included John M. Doe, Computer
Science, Spring 2015, sample committee/co-advisor, abstract, copyright,
acknowledgements, dedication to a dog, and lists of tables/figures. These
examples were identified as template material and preserved in full in docs/.

Class-controlled formatting remains untouched:
  - Letter paper, 12pt, openany, twoside; Times text with T1 encoding.
  - 1-inch nominal margins, heightrounded geometry, 0.5-inch footskip.
  - Double spacing, plain page style, bottom page numbers, heading sizes and
    spacing, float/caption settings, footnotes, indentation, and widow/orphan
    penalties.
  - Title and abstract layouts, committee listing, TOC layout, copyright-page
    numbering suppression, Roman preliminary pages, Arabic main pages.
  - Bibliography formatting and appendix numbering after the bibliography.

The example manual says 'oneside', but the actual class forces 'twoside'.
The implementation was preserved rather than changing it to match the stale
example prose. The class has no separate signature-page command: committee
names are printed on the title page. It defines no mathematical theorems.
This setup preserves the supplied class; it does not certify that this 2020
class meets any university requirements introduced since it was published.

Metadata declarations are \title, \author, \email, \department, \semester,
\advisor, repeated \committee, \mycopyright, \abstract, and (when wanted)
\acknowledgements. \coadvisor is optional. Defaults are examples, so omission
alone can print incorrect sample metadata. Email is stored but not visibly
printed by the class. There are no separate university or degree-name setters;
the university is literal class text and the degree follows the class option.

Matter order is preserved:
  \frontmatter, \maketitle, \makemycopyright, \makeabstract, TOC;
  \mainmatter, chapter includes;
  \backmatter, BibTeX bibliography;
  \appendix, then any future appendix chapters.
Optional acknowledgements, dedication, and empty figure/table lists are not
printed. Their insertion points and the complete original examples remain
available. No fictional acknowledgements or dedication were retained.

Baseline compilation before edits
---------------------------------

The untouched main was compiled with pdfLaTeX in nonstop mode, followed by
BibTeX. It did NOT build successfully: pdfLaTeX exited 1 and BibTeX exited 2.
A 22-page PDF was emitted despite errors; an emitted PDF was not treated as
proof of successful compilation. The previously supplied PDF also contained
unresolved references and stale master's-template metadata.

Pre-existing problems:
  - Missing pso and banana sample graphics (five failed includegraphics calls).
  - Float(s) lost and 'This may be a LaTeX bug' in the landscape example.
  - Overfull sample-figure text boxes up to 28.32971pt.
  - Undefined table reference and five undefined sample citations.
  - \bibliography{sample}, although no sample.bib exists; no usable td.bbl.

The active document now uses bibiTD and omits the archived sample demonstrations.
No class patch was needed. The original source and its instructional comments
remain byte-for-byte in docs/template-example.tex.txt.

Mathematical header audit
-------------------------

The root headerVarillyDiff.tex and
502/CIMPA_ECCO_2026/BiggerSymFncsECCO/headerVarillyDiff.tex are byte-identical
(1,231 lines each). The root version is canonical.

The audit inspected definitions and body usage in all 16 non-TD LaTeX files
under DoctoralNotebook, distinguishing six sources in Papers/ from ten other
notebook/workspace files. Counts below exclude Papers and strip comments;
duplicated drafts and text inside \iffalse remain, so these are occurrence
counts rather than counts of independent mathematical uses. All 72 imported
commands have use evidence outside Papers/.

Key body-use evidence:
  \bA        33 occurrences
  \bC       176 occurrences
  \bP       544 occurrences
  \cM        48 occurrences
  \ov       473 occurrences
  \la       526 occurrences
  \ps       474 occurrences
  \cO       508 occurrences
  \glu      259 occurrences
  \La       203 occurrences
  \Spec     137 occurrences
  \DefB      89 occurrences
  \ev        52 occurrences
  \Aut       42 occurrences
  \Ob        32 occurrences
  \vir       14 occurrences


The main notebook DoctoralNotebook/DN/DoctoralNotebookJIRR.tex supplies the
principal evidence: stable and pseudostable curves, Hodge integrals, stable
maps, virtual classes, localization, and the candidate stack of pseudostable
maps near line 3422. Its Budapest abstract near lines 3730-3749 describes the
construction as ongoing work, which justifies the tentative third title.

Imported groups (all defined explicitly in VarillyThesisMacros.tex):
  Blackboard: bA bC bE bL bP bQ bR bZ.
  Calligraphic/categories: cC cF cG cI cL cM cO cP cT cU gS cat Sch Curves.
  Greek: al bt Dl dl Ga ga kp La la om sg vf ze.
  Symbols: ov un x ox less hookto onto isom To del.
  Delimiters/text: set Set bonj genr word defeq.
  Operators: Aut codim DefB End ev Ext glu Hom id Isom Mor Ob PGL Pic pr
             Spec Stab Sym tot.
  Superscripts: ps vir.

Dependencies are only amsmath[cmex10], amsthm, amssymb (the latter loads
amsfonts). These were already present in the template and now have a single
home in the reduced header. Familiar environment names Th, Prop, Lem, Cor,
Def, Ex, Rmk have English titles and share a counter reset by chapter. This
avoids numbering such as 1.0.1 before a first section. The class's existing
equation numbering is unchanged.

Classifications and deliberate exclusions:
  1. Safe mathematical shorthand: the evidence-based subset listed above.
  2. Operators/theorems: relevant geometry operators and seven standard
     amsthm environments; no class command conflicts were found.
  3. Required packages: only the three AMS packages above. The main retains
     the template's graphicx, subfig[caption=false], booktabs, url, and cite.
  4. Deferred conveniences: quot, half, bold symbols/bm, matrix shortcuts,
     extra alphabets, starred theorems, and general-purpose diagram packages.
     Add a narrowly justified dependency when actual chapter content needs it.
  5. Presentation material: omit course/title metadata, worksheets, MATLAB
     listings, algorithm packages, index helpers, proof/answer boxes, Spanish
     language/hyphenation setup, colors, dark themes, and decorative TikZ styles.
     Notebook-local Mickey/Submickey, fundclass, and AI-summary helpers were
     also omitted; their display conventions are not needed by the baseline.
  6. Conflicting formatting: omit geometry/savetrees, memoir-only layout
     commands, fancyhdr, section/title redefinitions, TOC settings, paragraph
     indentation and line spread. Omit newpxtext, mathalfa, and mathabx, so
     the existing class/math fonts and standard symbols remain in control.
     The notation uses ordinary \mathcal, not the header's Euler substitution.
     Do not combine the header's subcaption with the template's subfig.
     Do not import the header's hyperref metadata, bibliography/indices, or
     packages with unrelated global behavior. The template leaves hyperref
     optional; tocloft supplies its needed no-op \phantomsection without it.
  7. Dangerous/dependent definitions: omit redefinitions of \l, \., \:,
     \mid, \div, \vec, \leq, \geq, and package-dependent \2/\third.
     Preserve LaTeX accent commands for author names and bibliography text.
     The old \ev uses renewcommand after physics; without physics it would
     fail, so the thesis defines the evaluation operator directly.
     drawcirculararc reuses cA/cB as scratch variables; late internal-@ helper
     definitions and custom theorem-ending overrides are not imported.

The root header has no need to be loaded by this project. No geometry, font,
page-style, title/section formatting, color, or bibliography package is placed
in VarillyThesisMacros.tex. The original \cref convenience remains 'Chapter',
so loading cleveref later would require resolving that name explicitly.

Bibliography consolidation
--------------------------

Before: bibiTD.bib had 5 entries; bibiDoctoralNotebook.bib had 30.
Immediately after the initial merge: bibiTD.bib had 35 unique entries.
The five template samples were subsequently removed by the author; the current
file has 30 entries. This follow-up sorting/build change leaves that file
byte-for-byte unchanged. No exact duplicate keys, conflicting
keys, or reliably identifiable duplicate works were found by normalized title,
DOI, arXiv identifier, or ISBN. No key was removed or renamed.

Repairs in the destination only:
  VakilRisingSea: two title fields became one title,
    'The Rising Sea: Foundations of Algebraic Geometry'.
    Both original title strings are retained in this combined value.
  JeongseokVirtClasses: @inbook became @incollection, since its existing fields
    describe an authored chapter in an edited book. All fields were retained.

The source's duplicate-title warning and author/editor warning are resolved.
Accents, capitalization braces, URLs, DOI/arXiv fields, abstracts, and other
metadata were preserved in the initial merge. The current 30 entries pass
BibTeX with plain without warnings. The standard plain style does not display
every stored DOI, URL, or arXiv field;
those fields remain in the database for any later citation-style decision.

The temporary \nocite{*} in td.tex renders all 30 current references because the
chapter stubs contain no citations. This is a setup inventory, not a claim that
all sources belong in the finished dissertation. Remove that line when actual
chapter citations exist; BibTeX then prints the cited references only.
The plain style controls alphabetical order, independent of database order,
while keeping the numeric citation format. The supplied template permits other
bibliography styles. The removed sample keys were forney2011classify, cebl,
forney2011thesis, forney2015echostate, haykin2009neural.

Unverified metadata issues preserved rather than silently guessed:
  - VFCWorkMath spells an author 'Battistell, Luca'.
  - KapranovPaper contains \overline{M}\_{0,n} (escaped underscore in math).
  - VakilRisingSea has a literal 'doi:' prefix in the DOI field.
  - The original (now removed) sample cebl spelled 'Compter Science'.
These do not cause compilation errors; they need bibliographic checking before
submission. This merge is a syntax/deduplication audit, not a literature audit.

Metadata and remaining author input
-----------------------------------

Filled:
  Author: José Ignacio Rojas Rojas.
    Evidence: personal signed closing in 340I2026/S/S340I2026.tex:241;
    the instructor there and the canonical header consistently use Ignacio
    Rojas. The expanded signature supplies the full name without invention.
  Email: ir@colostate.edu, from the same syllabus's instructor listing.
  Department: Department of Mathematics.
    Evidence: DoctoralNotebook/TM/tm_cover.tex:17 also names Colorado State
    University and Fort Collins, supporting the supplied class's institution.
  Advisor: Renzo Cavalieri.
  Committee: Mark Shoemaker, Maria Gillespie, Joshua Berger.
    These four names follow the user's explicit information and use \advisor
    plus three \committee declarations. No co-advisor was invented.

Still needed: approved thesis title, actual abstract, defense/submission term
and year, copyright statement/license and year, and any desired dedication or
acknowledgements. Confirm the exact degree wording before submission; Doctor
of Philosophy follows the original [doctor] selection and doctoral context.
No date or degree was copied from the earlier master's thesis.
Visible TODO text is deliberately incomplete content, not a finished abstract,
title, or copyright declaration. The third chapter title is tentative.

Figures and other assets
------------------------

Created only figures/ (with .gitkeep); no tikz/ or tables/ subdirectory yet.
The main uses \graphicspath{{figures/}} and extensions .pdf/.png/.jpg/.jpeg.
graphicx is actually loaded once, including its indirect load by pdflscape in
the class. No absolute graphics paths or draft-mode missing-image substitutes
are used. No graphical asset is required by the chapter stubs, and no existing
figure was copied or moved. There are consequently no copied-asset provenance
records or current missing/uncited graphical dependencies.

Future filenames: descriptive lowercase with hyphens, no spaces, special
characters, or unnecessary chapter/figure numbers; e.g. localization-graph.pdf.
Prefer PDF for vectors, PNG for raster diagrams/screenshots, JPEG for photos,
and TikZ source for native LaTeX diagrams. Keep table code in its chapter until
large reusable tables/data justify a tables/ directory. Never store build
artifacts in figures/. Figure captions go below figures; table captions go
above tables, per the preserved template instructions.

Possible later sources (not imported):
  DoctoralNotebook/figs/FigsDNnotability1.pdf through FigsDNnotability5.pdf:
    curves, dual graphs, strata, Hodge/psi classes, and intersection drawings;
    notebook uses cropped regions of larger sheets.
  DoctoralNotebook/figs/figLineBundleDefn.pdf
  DoctoralNotebook/figs/GluingMapActionEx.jpeg
  DoctoralNotebook/figs/DiagramGenus1TotChernG1PB.jpeg
  DoctoralNotebook/TM/Figs/fig1Intro.pdf and fig2Intro.pdf
  DoctoralNotebook/TM/Figs/AutC.pdf, DefC.pdf, normalBdlPb.pdf, tangentBdlPb.pdf
  DoctoralNotebook/TM/Figs/standaloneFigMaker.tex:
    source exists, but currently renders only one fragment; not demonstrated
    to reproduce all the PDFs above.

These assets have no confirmed authorship/licensing record from the inspected
files. Presence in the notebook and Quartz/pdfTeX producer metadata alone do
not establish ownership. On future import, record the original path, creator,
source file, cropping, and any citation/license. Preserve originals. Avoid the
old TM source's 'figs/' versus actual 'Figs/' case mismatch on portable builds.

Verification
------------

The restructured main compiles from TD with pdfLaTeX (TeX Live 2025), BibTeX
0.99d, and latexmk's settling passes. The current result is a 10-page PDF:
title, copyright placeholder, abstract placeholder, TOC, three chapter stubs,
and the complete provisional bibliography. All 30 current entries are recognized.
All three chapter previews also build successfully under build/. BUILDING.md
records the current commands and snapshot limitations. The source bibliography,
class, mathematical header, and chapter source files are unchanged by this
sorting/build follow-up.
Final checks cover duplicate keys/labels, included chapters, actual graphics
package load count, paths, layout preservation, and the final LaTeX/BibTeX logs.

A separate temporary, self-contained copy exercises all 72 mathematical
commands, all seven theorem environments and proof, representative moduli/Hodge/
virtual-class notation, chapter and theorem references, a real citation, and
an actual temporary PDF loaded by its basename from figures/. Those check
expressions and temporary graphical assets are not inserted into thesis prose.

Final verification result: PASS. No final LaTeX or BibTeX warnings, errors,
undefined commands/citations/references, duplicate labels/keys, missing files,
or overfull/underfull boxes. The temporary comprehensive typesetting build
also has no warnings. A negative test confirms a missing figure causes an
explicit compilation failure. Original and final geometry dimensions match.
The class and archived original source match their pre-edit copies exactly.
Visual review of the title page, TOC, tentative third chapter, and bibliography
confirmed the expected layout and readable text, including the accented name.
