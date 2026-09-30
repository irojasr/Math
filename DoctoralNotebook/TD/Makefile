LATEXMK ?= latexmk

.DEFAULT_GOAL := thesis
.PHONY: thesis introduction background chapter1 clean
# Serialize shared full-build bootstrapping, including when make -j is used.
.NOTPARALLEL:

thesis:
	$(LATEXMK) td.tex

# Seed includeonly's labels and counter checkpoints from the last full build.
# Separate auxiliary directories keep previews from overwriting that snapshot.
introduction background chapter1:
	@for aux in build/td.aux build/chapters/introduction.aux \
		build/chapters/background.aux build/chapters/chapter1.aux; do \
		if test ! -s "$$aux"; then $(MAKE) thesis || exit $$?; break; fi; \
	done
	@mkdir -p build/previews/$@/chapters
	@cp build/td.aux build/previews/$@/$@-preview.aux
	@cp build/chapters/*.aux build/previews/$@/chapters/
	$(LATEXMK) -g -outdir=build/previews/$@ -out2dir=build \
		-jobname=$@-preview \
		-usepretex='\def\ThesisPreviewChapter{chapters/$@}' td.tex

clean:
	rm -rf -- build
