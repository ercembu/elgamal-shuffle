#!/usr/bin/bash

cd ${PWD}/thesis

pdflatex thesis
bibtex thesis
pdflatex thesis
