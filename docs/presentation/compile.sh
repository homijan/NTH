#! /bin/sh
rm -r build
mkdir build && cd build

cp ../MultiMat2017.tex document.tex
#cp ../CELIA.tex document.tex

pdflatex document.tex
pdflatex document.tex

mv document.pdf ../MultiMat2017.pdf
#mv document.pdf ../CELIA.pdf
cd ..
rm -r build
