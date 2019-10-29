# Adapted makefiles for feap (8.4.1d/e)
# Version 1.1

# This makefile allows parallel building ( make -j #, # = number of your processor cores)

# In this file you find the standard object files
include objects.mk

# Which compiler, 64 or 32 bit?
COMPILER=MY
#COMPILER=64
#COMPILER=32


# This is for you 
MYOBJS= 
#./myObjs/user_material.o ./myObjs/CalcInverse.o
#./damageRoutines/elmt05.o   ./damageRoutines/e03help.o   ./damageRoutines/elmthelp.o   ./damageRoutines/elmtmat.o   ./damageRoutines/gaussmod.o   ./damageRoutines/gaussscha_1.o


ifeq ($(COMPILER), MY)
# Here you can choose your own compiler and flags
FC=ifort
CC=gcc
FFLAGS=-O2 -g
CFLAGS=-pipe -O2 -g
# for 32 bit
#INCLUDE=-I../../include -I../../include/integer4 -I../../include/openmpi
# for 64 bit
INCLUDE=-I../../include -I../../include/integer8 #-I../../include/openmpi
else
# These are standard (GNU) compilers, do not modify
ifeq ($(COMPILER), 64)
FC=gfortran
CC=gcc
FFLAGS=-pipe -O2 -g -m64
CFLAGS=-pipe -O2 -g -m64
INCLUDE=-I../../include -I../../include/integer8 #-I../../include/openmpi
else
FC=gfortran
CC=gcc
FFLAGS=-pipe -O2 -g -m32
CFLAGS=-pipe -O2 -g -m32
INCLUDE=-I../../include -I../../include/integer4 
#-I../../include/openmpi
endif
endif

# without MKL
# LFLAGS=-lm -lX11 -llapack
# with MKL
LFLAGS=-lm -lX11 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread



# Do not touch this! Here be dragons!

OBJECTS=$(STDOBJECTS) $(MYOBJS)

.SUFFIXES : .o .c .f .f90

feap: $(OBJECTS)
	$(FC) $(OBJECTS) $(LFLAGS) -o feap
	@echo
	@echo Hooray,,1! Feap was built successfully
	@echo

feaplib: $(OBJECTS)
	ar rv libfeap.a $(OBJECTS)
	@echo
	@echo Hooray,,1! Feap library was built successfully
	@echo

clean:
	-rm -f $(OBJECTS)

allclean: clean
	-rm -f feap libfeap.a

veryclean: 
	find ../../ -name '*.o' -exec rm {} \;
	-rm -f feap libfeap.a

mrproper:
	find ../../ -name '*.o' -exec rm {} \;
	find ../../ -name 'built_objs' -exec rm {} \;
	-rm -f feap libfeap.a ../../Feap8_4.a ../../main/feap

help:
	@echo Usage:
	@echo Call make with -j\# \(\# being the number of your processor cores, i.e. a parallel build\)
	@echo 
	@echo Targets \(make TARGET\)
	@echo "feap         - Default target, builds feap"
	@echo "feaplib      - builds libfeap.a"
	@echo "clean        - deletes all object files (all objects mentioned in "
	@echo "               Makefile and objects.mk, others are ignored)"
	@echo "allclean     - performs clean and deletes feap"
	@echo "veryclean    - delete every *.o file via find. USE WITH CAUTION!"
	@echo 
	@echo Don\'t Panic

# These are internal targets and should never becalled directly

.f.o:
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $*.o $*.f

.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $*.o $*.c

.f90.o:
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $*.o $*.f90
