include ../makefile.in

OBJECTS = stripm.c

stripm:
	$(CC) -o stripm $(OBJECTS)

.c.o:
	$(CC) -c $(CCOPTFLAG) $*.c -o $*.o

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.c,v

rcs:
	mkdir RCS
	ci -t-"" *.c
	ci -t-"" makefile

