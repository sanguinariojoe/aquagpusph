#!/bin/bash
# To generate the PDF you need to install:
# 1	gnuplot
# 2	python
# 3	LaTeX

# Plot some images
gnuplot images/lateral_water_1x_deleffe/lateral_water_1x.gnuplot
python images/wendland/wendland2D.py

# Compile
pdflatex -interaction=nonstopmode main.tex
bibtex main.aux
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex

# Rename output file
mv main.tex aquagpusph-usermanual.pdf