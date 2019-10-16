include ../../makefile.in

clean: 
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile

