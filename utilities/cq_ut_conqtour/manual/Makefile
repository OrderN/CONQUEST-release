MANUAL=ConqtourManual
TMP=$(MANUAL).aux\
    $(MANUAL).idx\
    $(MANUAL).ilg\
    $(MANUAL).ind\
    $(MANUAL).log\
    $(MANUAL).out\
    $(MANUAL).toc\
    $(MANUAL).blg\
    $(MANUAL).bbl

default: $(MANUAL).pdf

$(MANUAL).pdf: $(MANUAL).tex
	pdflatex $(MANUAL).tex
#	bibtex $(MANUAL)
	makeindex $(MANUAL)
	pdflatex $(MANUAL).tex
	pdflatex $(MANUAL).tex

clean:
	rm -f $(TMP)
	
deepclean:
	rm -f $(TMP) $(MANUAL).pdf