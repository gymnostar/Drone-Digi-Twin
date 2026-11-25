#!/bin/sh
# Thursday 26 Nov deck: numbers from the JSON results -> numbers.tex -> slides.pdf
#   sh DLR/slides/2026-11-26/build.sh
cd "$(dirname "$0")" || exit 1
REPO=$(cd ../../.. && pwd)
PATH="/Library/TeX/texbin:$PATH"; export PATH
"$REPO/.venv/bin/python" make_numbers.py || exit 1
latexmk -pdf -interaction=nonstopmode -halt-on-error -silent slides.tex > build.log 2>&1 \
  && { latexmk -c slides.tex >/dev/null 2>&1; rm -f slides.nav slides.snm slides.vrb build.log; echo "$(pwd)/slides.pdf"; } \
  || { echo "LaTeX FAILED, see build.log"; tail -25 build.log; exit 1; }
