FIGFILES = $(wildcard *.fig)
AGRFILES = $(wildcard *.agr)
TARGET   = $(FIGFILES:%.fig=%.pdf) $(AGRFILES:%.agr=%.pdf)

all: $(TARGET)

veryall: all

.PHONY: clean veryclean

%.eps:	%.fig
	fig2dev -L eps $< $@

%.pdf:	%.fig
	fig2dev -L pdf $< $@

%.png:	%.fig
	fig2dev -L png $< $@

%.pdf:	%.dot
	dot -Tpdf -o $@ $<

%.eps:  %.agr
	gracebat -hdevice EPS $*.agr

%.pdf:  %.agr
	gracebat -hdevice EPS $*.agr && \
	epstopdf $*.eps && \
	rm -f $*.eps

clean:
	-rm -f *.bak

veryclean: clean
	-rm -f $(TARGET)
