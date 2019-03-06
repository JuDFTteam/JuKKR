PROJECT=report
TEX=pdflatex
BIBTEX=bibtex
BUILDTEX=$(TEX) $(PROJECT).tex

all:
	    $(BUILDTEX)
		    $(BUILDTEX)
		        $(BUILDTEX)

open:
	    make
	    open $(PROJECT).pdf

clean-all:
	    rm -f *.dvi *.log *.bak *.aux *.bbl *.blg *.idx *.ps *.eps *.pdf *.toc *.out *blx.bib *.lof *run.xml *~

clean:
	    rm -f *.log *.bak *.aux *.bbl *.blg *.idx *.toc *.out *blx.bib *.lof *run.xml *~

