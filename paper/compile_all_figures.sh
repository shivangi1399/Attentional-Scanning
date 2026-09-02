#!/bin/bash
# =====================================================================
# Concatenate the paper figures into writing/figures/all_figures.pdf.
#
# Run this after figures.m, whenever any figure changes.
#
# ORDER is figures.m's block order, with 6b following 6. It is recorded
# here because it is recorded nowhere else -- the combined PDF used to be
# assembled by hand, so the sequence lived only in whoever last built it.
# If you add a figure to figures.m, add it here in the same position.
#
# pdfunite, not gs or a LaTeX wrapper: it preserves each figure's own page
# box. The figures range from 15x15 cm to 19x27.5 cm and are deliberately
# not scaled onto a common sheet, so a tool that imposes one page size
# would silently crop or letterbox them.
# =====================================================================
set -o pipefail

BASE=${BASE:-/mnt/hpc/projects/MWSampling/4Shivangi}
P=$BASE/writing/figures/plots
OUT=$BASE/writing/figures/all_figures.pdf

PAGES=(
    "$P/fig03_erp_grandavg.pdf"
    "$P/fig04_tfr.pdf"
    "$P/fig05_phase_measures.pdf"
    "$P/fig06_hypotheses.pdf"
    "$P/fig06b_hypothesis_differences.pdf"
    "$P/fig07a_phase_progression.pdf"
    "$P/fig07b_pgd.pdf"
    "$P/fig07c_derotation.pdf"
    "$P/fig07d_per_position_hermes.pdf"
    "$P/fig07d_per_position_klecks.pdf"
)

# Refuse to build a short PDF from a missing figure: pdfunite would happily
# skip it and the gap is easy to miss in a 10-page document.
missing=0
for f in "${PAGES[@]}"; do
    [ -f "$f" ] || { echo "MISSING: $f"; missing=1; }
done
[ "$missing" -eq 0 ] || { echo "aborted -- run figures.m first"; exit 1; }

pdfunite "${PAGES[@]}" "$OUT" || { echo "pdfunite failed"; exit 1; }

echo "wrote $OUT"
pdfinfo "$OUT" | grep -E '^Pages:'
echo "expected ${#PAGES[@]} pages"
