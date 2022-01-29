main.pdf: main.tex imports/titlepage.tex pre_clean prepare_bibliography
	pdflatex --shell-escape main.tex
	bibtex main
	pdflatex --shell-escape main.tex
	pdflatex --shell-escape main.tex
	cp main.pdf ../thesis.pdf
	$(MAKE) pre_clean

prepare_bibliography:
	python bib/prepare_bib.py

pre_clean:
	rm -f *.aux *.blg *.log *.bbl *.toc *.out *.bcf *.xml
	rm -f imports/*.aux
	rm -rf _minted-main

clean: pre_clean
	rm -f *.pdf
