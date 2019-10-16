include ../makefile.in

SOURCES = $(FSOURCE)*.$(FEXT)

archive: $(SOURCES)
	$(AR) $(ARFEAP) $%
	@@echo "--->WINDOWS objects now up-to-date<---"
	rm -f $% *.o

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v

install: $(SOURCES)
	@@echo "--->WINDOWS objects not installed in archive<---"

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile

.f.o:
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f -o $*.o

.c.o:
	$(CC) -c $(CCOPTFLAG) *.c -o $*.o

