include ../../makefile.in

SOURCES = *.f

OBJECTS = *.o

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v

install: $(SOURCES)
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f
	ar -rv arpacklib.a *.o
	ranlib arpacklib.a
	@@echo "--->ARPACK objects now installed in archive<---"
	rm -f $% *.o

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile

.f.o:
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f -o $*.o

.c.o:
	$(CC) -c -I$(CINCLUDE) $(CCOPTFLAG) *.c -o $*.o

