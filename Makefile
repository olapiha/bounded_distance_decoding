thesis:
	pdflatex $(file)
	bibtex $(file)
	pdflatex $(file)
	pdflatex $(file)
	make clean

clean:
	@rm -f *.aux
	@rm -f *.bcf
	@rm -f *.idx
	@rm -f *.log
	@rm -f *.out
	@rm -f *.bbl
	@rm -f *.blg
	@rm -f *.xml
	@rm -f *.toc
	@rm -f *.fls
	@rm -f *.fdb_latexmk
	@rm -f *.synctex.gz
	@rm -f __latexindent*
