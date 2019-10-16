include ../makefile.in

OBJECTS = feap84.o

SOURCES = $(FSOURCE)feap84.$(FEXT)

feap: $(OBJECTS) $(SOURCES) $(ARFEAP)
	ranlib $(ARFEAP)
	$(FF) -o feap $(OBJECTS) $(ARCHIVELIB) $(ARPACKLIB) \
	$(BLASLIB) $(LAPACKLIB) $(ARFEAP) $(LDOPTIONS)

.f.o:
	$(FF) -c $(FFOPTFLAG) -I$(FINCLUDE) $*.f -o $*.o

.c.o:
	$(CC) -c $(CCOPTFLAG)  $*.c -o $*.o

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile

