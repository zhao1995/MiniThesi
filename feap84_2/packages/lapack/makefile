
include ../../makefile.in

SOURCES = *.f

clean: 
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v

install: $(SOURCES)
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f
	ar -rv lapacklib.a *.o
	ranlib lapacklib.a
	@@echo "--->LAPACK objects now installed in archive<---"
	rm -f $% *.o

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile

.f.o:
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f -o $*.o

.c.o:
	$(CC) -c $(CCOPTFLAG) *.c -o $*.o

