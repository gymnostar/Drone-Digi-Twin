#!/bin/sh
# One command for the meeting material: tables + pictures + PDF(s). Works with partial results.
#   sh DLR/meshStudy/report/build.sh                    (from anywhere)
#   SKIP_FIGURES=1 sh .../build.sh                      skip the ParaView pictures (fast; keeps the old ones)
#   sh .../build.sh --include-running                   also show probe values of unfinished runs (marked)
# Steps: make_table.py (md/html + 2 plots) -> make_figures.py (pictures, pvpython) -> make_tex.py
#        (tables.tex, summary.tex, fig/) -> pdflatex via latexmk -> mesh_study.pdf (+ slides.pdf)
cd "$(dirname "$0")" || exit 1
REPORT=$(pwd)
REPO=$(cd ../../.. && pwd)
PY="$REPO/.venv/bin/python"
PATH="/Library/TeX/texbin:$PATH"; export PATH
[ -x "$PY" ] || { echo "python venv not found at $PY"; exit 1; }

echo "== 1/4 make_table.py"
"$PY" ../make_table.py > make_table.out 2>&1 \
  && echo "   ok (output in report/make_table.out)" \
  || echo "   make_table.py FAILED (see report/make_table.out); continuing with its old plots"

if [ "${SKIP_FIGURES:-0}" = "1" ]; then
  echo "== 2/4 make_figures.py skipped (SKIP_FIGURES=1)"
else
  echo "== 2/4 make_figures.py (ParaView cuts, one case at a time)"
  "$PY" ../make_figures.py || echo "   make_figures.py FAILED; continuing with the pictures that exist"
fi

echo "== 3/4 make_tex.py"
"$PY" make_tex.py "$@" || { echo "make_tex.py FAILED"; exit 1; }

echo "== 4/4 LaTeX"
build_pdf() {
  [ -f "$1.tex" ] || return 0
  rc=1
  if command -v latexmk >/dev/null 2>&1; then
    latexmk -pdf -interaction=nonstopmode -halt-on-error -silent "$1.tex" > "$1.build.log" 2>&1; rc=$?
  elif command -v pdflatex >/dev/null 2>&1; then
    pdflatex -interaction=nonstopmode -halt-on-error "$1.tex" > "$1.build.log" 2>&1 \
      && pdflatex -interaction=nonstopmode -halt-on-error "$1.tex" >> "$1.build.log" 2>&1; rc=$?
  else
    echo "   no LaTeX engine found: zip $REPORT and compile $1.tex on Overleaf (pdfLaTeX)"; return 1
  fi
  if [ $rc -eq 0 ] && [ -f "$1.pdf" ]; then
    echo "   $REPORT/$1.pdf"
    command -v latexmk >/dev/null 2>&1 && latexmk -c "$1.tex" >/dev/null 2>&1
    rm -f "$1.nav" "$1.snm" "$1.vrb" "$1.out" "$1.toc"
    rm -f "$1.build.log"
  else
    echo "   $1.tex FAILED to compile, see $REPORT/$1.build.log"; tail -20 "$1.build.log"; return 1
  fi
}
build_pdf mesh_study
build_pdf slides
echo "done."
