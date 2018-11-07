# LaTeX Makefile
FILE=manuscript
all: $(FILE).pdf \

docx: $(FILE).docx

# clean up processing files
.PHONY: clean
clean:
	\rm *.aux *.blg *.out *.bbl *.log *.spl *.log

# PDF export
$(FILE).pdf: $(FILE).tex \
	header.tex
	nihunsrt.bst \
	bibliography.bib
	pdflatex $(FILE)
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)
#	pdftk $(FILE).pdf update_info report.txt output $(FILE)_copy.pdf
#	mv $(FILE)_copy.pdf $(FILE).pdf

# DOCX export
$(FILE).docx: $(FILE).tex
	cp $(FILE).tex $(FILE)_docx_reformat.tex
	sed -i 's/pdf/png/g' $(FILE)_docx_reformat.tex # use png versions of figures
	sed -i 's/\\(\\kappa\\)/κ/g' $(FILE)_docx_reformat.tex
	sed -i 's/\\(\\beta\\)/β/g' $(FILE)_docx_reformat.tex
	sed -i 's/\\(\\alpha\\)/α/g' $(FILE)_docx_reformat.tex
	sed -i 's/\\(\\mu\\)/μ/g' $(FILE)_docx_reformat.tex
	sed -i 's/\\(\\gamma\\)/γ/g' $(FILE)_docx_reformat.tex
	sed -i 's/{\"i}/ï/g' $(FILE)_docx_reformat.tex
	sed -i 's/\\pm/±/g' $(FILE)_docx_reformat.tex
	sed -i 's/\num{//g' $(FILE)_docx_reformat.tex
	sed -i 's/\\(_{\\text{2}}\\)/₂/g' $(FILE)_docx_reformat.tex
	sed -i 's/\\(^{\\text{2}}\\)/²/g' $(FILE)_docx_reformat.tex
	pandoc --bibliography=bibliography.bib --filter pandoc-citeproc --reference-docx=template.docx --csl=science.csl --latex-engine=pdflatex --write=docx $(FILE).tex -o $(FILE).docx
