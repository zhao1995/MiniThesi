include ../makefile.in

SOURCES = $(wildcard $(FSOURCE)*.$(FEXT)) $(wildcard $(CSOURCE)*.$(CEXT))
OBJSF   = $(SOURCES:.f=.o)
OBJS    = $(OBJSF:.c=.o)
FL      = built_objs

.SILENT:

archive:
	@rm -f $(FL); \
	$(MAKE) $(ARFEAP); \
	MKRET=$$?; \
	if [ -f $(FL) ]; then \
	  $(AR) $(ARFEAP) `cat $(FL)`; \
	  rm `cat $(FL)`; \
	  echo "--->`pwd` objects now up-to-date<---"; \
	fi; \
	if [ $$MKRET -ne 0 ]; then false; fi;

$(ARFEAP): $(ARFEAP)($(OBJS))

clean: 
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile

.f.a:
	@if $(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) $*.f -o $*.o; then echo $*.o >> $(FL); else false; fi;
	@echo "Building $*.o"

.c.a:
	if $(CC) -c $(CCOPTFLAG) $*.c -o $*.o; then echo $*.o >> $(FL); else false; fi;
