#! /bin/sh
rm -r build
mkdir build && cd build

#cp ../NLTWDFEOS.tex document.tex
cp ../MultiMat2017.tex document.tex
#cp ../LANL2017.tex document.tex

pdflatex document.tex
pdflatex document.tex

mv document.pdf ../MultiMat2017.pdf
cd ..
rm -r build
