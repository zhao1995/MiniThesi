#  Compile the parallel parmetis partioner with 'make partition'
#  Run the partioner with 'make partrun'
#  See further notes below!

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

partition: partition.o sparse-graph.o
	-${CLINKER} -o partition partition.o sparse-graph.o \
		${PETSC_LIB}

# Set as need or just type in at the command line
NPROC=2
NPART=2
INAME=Iblockflat

#
# -n number controls the number of processors you will be using to partition
# your problem.   The number before the inputfile name controls the number of
# processors FOR WHICH you are partitioning the mesh
#
# The output will be graph.block (in the case below, strips off the leading I)
# Note that Iblock (in the case below) must be a flat file
#
# Next, start feap in the parfeap directory and then produce the parallel input
# files as:
# 
# Macro> graph,file 
# Macro> outd
# 
# Note 'file' is the word "file" not "graph.block".
#

partrun:
	-@${MPIEXEC} -n $(NPROC) \
	$(FEAPHOME8_4)/parfeap/partition/partition $(NPART) $(INAME)


rcs:
	mkdir RCS
	ci -t-"" *.h
	ci -t-"" *.c
	ci -t-"" makefile
	ci -t-"" README

checkout:
	co -q RCS/*.h,v
	co -q RCS/*.c,v
	co -q RCS/README,v

