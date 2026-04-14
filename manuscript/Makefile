# ---------------------------------------------------------------------------
# Makefile for the preprint LaTeX build
# Target platform: macOS with TeX Live 2026 (pdflatex + bibtex + latexmk)
# ---------------------------------------------------------------------------

TEX       := preprint.tex
BASE      := preprint
PDFLATEX  := pdflatex -interaction=nonstopmode -halt-on-error
BIBTEX    := bibtex

AUX_EXTS  := aux log bbl blg out toc lof lot fls fdb_latexmk synctex.gz run.xml bcf

# Default target
.PHONY: all
all: pdf

# ---------------------------------------------------------------------------
# Primary build: pdflatex -> bibtex -> pdflatex -> pdflatex
# ---------------------------------------------------------------------------
.PHONY: pdf
pdf: $(BASE).pdf

$(BASE).pdf: $(TEX) references.bib
	$(PDFLATEX) $(TEX)
	$(BIBTEX)   $(BASE) || true
	$(PDFLATEX) $(TEX)
	$(PDFLATEX) $(TEX)

# ---------------------------------------------------------------------------
# bioRxiv-ready build: single spacing, no line numbers
# Passes \def\biorxiv{1} on the command line before \input-ing the main tex.
# Produces preprint_biorxiv.pdf.
# ---------------------------------------------------------------------------
.PHONY: biorxiv
biorxiv: preprint_biorxiv.pdf

preprint_biorxiv.pdf: $(TEX) references.bib
	$(PDFLATEX) -jobname=preprint_biorxiv "\def\biorxiv{1}\input{$(TEX)}"
	$(BIBTEX)   preprint_biorxiv || true
	$(PDFLATEX) -jobname=preprint_biorxiv "\def\biorxiv{1}\input{$(TEX)}"
	$(PDFLATEX) -jobname=preprint_biorxiv "\def\biorxiv{1}\input{$(TEX)}"

# ---------------------------------------------------------------------------
# Live preview via latexmk
# ---------------------------------------------------------------------------
.PHONY: watch
watch:
	latexmk -pdf -pvc -interaction=nonstopmode $(TEX)

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
.PHONY: clean
clean:
	@for ext in $(AUX_EXTS); do \
	    rm -f $(BASE).$$ext preprint_biorxiv.$$ext; \
	done
	@rm -f *.auxlock

.PHONY: distclean
distclean: clean
	@rm -f $(BASE).pdf preprint_biorxiv.pdf
