# Springer Nature LaTeX template setup

The official *Nature Communications* LaTeX template is **`sn-jnl.cls`**
with the `nature` option. It is **not** distributed via CTAN or TeX
Live; Springer Nature publishes it directly. Two distribution channels:

1. **Overleaf template (one-click clone):**
   <https://www.overleaf.com/latex/templates/springer-nature-latex-template/myxmhdsbzkyd>
2. **Direct download (zip):** Springer Nature LaTeX author support
   <https://www.springernature.com/gp/authors/campaigns/latex-author-support>
   (look for "LaTeX template for Nature journals" — yields a zip with
   `sn-jnl.cls`, `sn-nature.bst`, sample BibTeX, sample tex).

## Install steps (local TeX Live / TinyTeX)

```bash
# 1. Download the zip from the Springer support page (one-time, manual).
# 2. Unpack into the manuscript directory so the .cls is visible to pdflatex:
mkdir -p manuscript/templates/sn-jnl
unzip ~/Downloads/sn-latex.zip -d manuscript/templates/sn-jnl

# 3. Tell pdflatex to look there:
export TEXINPUTS=./templates/sn-jnl:$TEXINPUTS

# 4. Verify:
kpsewhich sn-jnl.cls
```

Alternative: copy `sn-jnl.cls`, `sn-nature.bst`, and any required font
files directly next to `preprint_ncomms.tex`.

## Switch the build to use it

Open `manuscript/preprint_ncomms.tex` and flip the toggle near the top:

```tex
\snjnltrue   % was \snjnlfalse (community fallback)
```

Then rebuild:

```bash
cd manuscript
pdflatex preprint_ncomms.tex
bibtex   preprint_ncomms
pdflatex preprint_ncomms.tex
pdflatex preprint_ncomms.tex
```

## Stand-in (works today, no extra install)

If you don't want to fetch the official template right now, the
fallback uses the community `nature.cls` (already installed via
TinyTeX). It produces a Nature-styled PDF that is close enough for
internal review and Yueh-Hua's read-through. Editorial submission
should use the official `sn-jnl.cls`.

## Body content — section split

The plan is to split the body of `preprint.tex` into
`manuscript/sections/{background,results,discussion,methods,conclusions}.tex`
so both builds (`preprint.tex` bioRxiv build and
`preprint_ncomms.tex` submission build) share the prose. This is
deliberately deferred until prose polish is settled — otherwise the
splits get re-churned on every editorial pass.

Recommended trigger for doing the split: when the word-count target
(≤5000 words main text + Methods) is hit and Yueh-Hua signs off on
the prose. Run:

```bash
# proposed helper script (not yet implemented)
python scripts/split_sections.py manuscript/preprint.tex manuscript/sections/
```

## What the submission portal will accept

Springer Nature's submission system **does** accept compiled PDFs, and
will re-typeset on acceptance — so the *.cls* used at submission is
primarily a courtesy that produces a near-final-look PDF for review.
It is not technically required to use `sn-jnl.cls` at first submission;
the community `nature.cls` produces an acceptable PDF. The official
class becomes mandatory only at the accepted-manuscript stage.
