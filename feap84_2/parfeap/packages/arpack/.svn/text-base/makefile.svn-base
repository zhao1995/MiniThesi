
# FEAP include locations and directives
  include ../../../makefile.in


# For ver 3.1-3.4 PETSc installations use:
  include ${PETSC_DIR}/conf/variables
  include ${PETSC_DIR}/conf/rules
  PETSC_FC_INCLUDES := $(PETSC_FC_INCLUDES) -I$(FINCLUDE) -I../.. # N.B. := avoids recursion


SOURCES = parfeaps.F parkv.F parmv.F paropk.F paropm.F pminvsqr.F \
	pddot.f pdgetv0.f pdnrm2.f pdsaitr.f pdsaup2.f pdsaupd.f \
	pdseupd.f umacr4.f

OBJECTS = parfeaps.o parkv.o parmv.o paropk.o paropm.o pminvsqr.o \
	pddot.o pdgetv0.o pdnrm2.o pdsaitr.o pdsaup2.o pdsaupd.o \
	pdseupd.o umacr4.o

#	$(FF) -c $(FFOPTFLAG) $(PETSC_FC_INCLUDES) *.f *.F
install: $(OBJECTS)
	ar -rv parpacklib.a *.o
	ranlib parpacklib.a
	@@echo "--->PARPACK objects now installed in archive<---"
	rm -f $% *.o

rcsclean: 
	echo $(FEAPHOME8_4)
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v
	co -q RCS/*.F,v
	co -q RCS/*.h,v

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" *.F
	ci -t-"" *.h
	ci -t-"" makefile

.f.o:
	$(FF) -c $(FFOPTFLAG) -I$(FINCLUDE) $*.f -o $*.o

.c.o:
	$(CC) -c $(PETSC_CC_INCLUDES) $(CCOPTFLAG) $*.c -o $*.o

