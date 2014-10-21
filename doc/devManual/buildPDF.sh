#!/bin/bash
# To generate the PDF you need to install:
# 2	python
# 3	LaTeX

CURRENT_DIR=`pwd`

# Plot some images
cd ../userManual
python images/tree/tree.py > tree_py.tex

cd $CURRENT_DIR
cp ../userManual/tree_py.tex ./

# Compile
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex

# Rename output file
mv main.pdf aquagpusph-devmanual.pdf
