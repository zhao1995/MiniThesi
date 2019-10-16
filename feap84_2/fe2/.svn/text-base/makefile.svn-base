include ../makefile.in

OBJECTS = ../main/feap84.o \
	parsend.o \
	parstop.o \
	pform.o \
	pforma.o \
	pstart.o \
	rvemat.o \
	rvesr.o \
	srvemat.o \
	tinput.o \
	trvemat.o \
	umacr9.o \
	umacr10.o \
	umacr11.o \
	umesh3.o \
	umesh4.o \
	uprob.o \
	usetio.o

SOURCES = ../main/$(FSOURCE)feap84.$(FEXT) \
	parsend.f \
	parstop.f \
	pform.f \
	pforma.f \
	pstart.f \
	rvemat.f \
	rvesr.f \
	srvemat.f \
	tinput.f \
	trvemat.f \
	umacr9.f \
	umacr10.f \
	umacr11.f \
	umesh3.f \
	umesh4.f \
	uprob.f \
	usetio.f

install: $(OBJECTS) $(SOURCES) $(ARFEAP)
	ranlib $(ARFEAP)
	$(FFE2) -o feap $(OBJECTS) $(ARFEAP) $(LDOPTIONSFE2)


.f.o:
	$(FFE2) -c $(FFOPTFLAG) -I$(FINCLUDE) $*.f -o $*.o

.c.o:
	$(CCE2) -c $(CCOPTFLAG)  $*.c -o $*.o

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v
	co -q RCS/*.h,v

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" *.h
	ci -t-"" makefile

