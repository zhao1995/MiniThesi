include ../../makefile.in

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.h,v

rcs:
	mkdir RCS
	ci -t-"" *.h
	ci -t-"" makefile
