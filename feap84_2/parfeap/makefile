# Set number of processors
  NPROC     = 2

# Set path options for executable
 FEAPRUN   = $(FEAPHOME8_4)/parfeap/feap
#
# Use this variant of FEAPRUN when debugging with valgrind
#
# FEAPRUN   = valgrind -v --leak-check=full --show-reachable=yes $(FEAPHOME8_4)/parfeap/feap

# The following two library package can be omitted if you wish to omit the
# the parallel eigensolve capability
#
# arpack archive library
  ARCHIVELIB = $(FEAPHOME8_4)/packages/arpack/archive/archivelib.a
# parpack library
  PARPACKLIB = $(FEAPHOME8_4)/parfeap/packages/arpack/parpacklib.a

# Local lapack library (use if you do not have lapack already)
#  LAPACKLIB = $(FEAPHOME8_4)/packages/lapack/lapacklib.a

# Local blas library
#  BLASLIB   = $(FEAPHOME8_4)/packages/blas/blaslib.a    

# FEAP make directives
  include $(FEAPHOME8_4)/makefile.in

# PETSc include files
  include ${PETSC_DIR}/conf/variables
  include ${PETSC_DIR}/conf/rules

  PETSC_FC_INCLUDES := $(PETSC_FC_INCLUDES) -I$(FINCLUDE)    # N.B. := avoids recursion

OBJECTS = ../main/feap84.o \
	adomnam.o filargs.o fparsop.o parbmat.o pcompress.o \
	pdomain.o pfeapmi.o pfeapsr.o pform.o pnodal.o parstop.o parmacio.o \
	pmacr7.o unix/pmetis.o unix/smetis.o pmodify.o ppmodin.o pparlo.o \
	parform.o pndata.o prwext.o psetb.o psproja.o psprojb.o pstart.o \
	psubsp.o scalev.o uasble.o uasblem.o upremas.o uassdk.o uassdm.o \
	usolve.o upc.o pextnd.o tinput.o pprint.o uldout.o \
        constass.o crsurfinv.o cpoutm.o cpoutmc.o cpoutms.o

install: $(OBJECTS) $(ARFEAP)
	ranlib $(ARFEAP)
	-${FLINKER} ${FFLAGS} -o feap $(OBJECTS) \
	$(PARPACKLIB) \
	$(ARCHIVELIB)  $(ARFEAP) \
	${PETSC_FORTRAN_LIB} ${PETSC_LIB} $(LDOPTIONS)
	@@echo "--->PARFEAP:feap now built"
#	rm -f $% *.o

checkout:
	co -q RCS/*.h,v
	co -q RCS/*.f,v
	co -q RCS/*.F,v

rcs:
	mkdir RCS
	ci -t-"" *.h
	ci -t-"" *.f
	ci -t-"" *.F
	ci -t-"" makefile

rcsclean:
	rcsclean -q
	rm -f *.o

#------------------- TYPICAL JOB SUBMISSIONS-----------------------------------
# Standard mpirun submissions works with AIJ and BAIJ matrix formats
feaprun:
	-@${MPIEXEC} -n $(NPROC)  $(FEAPRUN) -ksp_type cg -ksp_monitor -pc_type jacobi  -memory_info -malloc_log -log_summary  -options_left

feaprun-gmres:
	-@${MPIEXEC} -np $(NPROC) $(FEAPRUN) -ksp_type gmres  -pc_type bjacobi -log_summary -on_error_attach_debugger

# Prometheus (N.B. Prometheus is being fased out and replaced by GAMG).  Requires OUTDomain,BAIJ,,<block_size> matrix format.
feaprun-mg:
	-@${MPIEXEC}  -n $(NPROC)  $(FEAPRUN) -ksp_type cg -ksp_monitor -log_summary -out_verbose 2 -pc_type prometheus -pc_mg_type multiplicative -pc_mg_cycles 1 -prometheus_preduce_base 500 -prometheus_top_grid_limit 2500 -aggmg_avoid_resmooth -aggmg_smooths 1  -prometheus_repartition -options_left

# Hypre-Boomer example.  Requires OUTDomain,AIJ,<1>,<block_size>
feaprun-hypre-boomer:
	-@${MPIEXEC} -n $(NPROC) $(FEAPRUN) -ksp_type cg -ksp_monitor -pc_type hypre -pc_hypre_type boomeramg -log_summary -pc_hypre_boomeramg_strong_threshold 0.25 -pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi -ksp_view -on_error_attach_debugger -options_left

# ML example.  Requires OUTDomain,AIJ,<1>,<block_size>
feaprun-ml:
	-@${MPIEXEC} -n $(NPROC) $(FEAPRUN) -ksp_type cg -ksp_monitor -pc_type ml  -mat_no_inode -log_summary -ksp_view_binary  -options_left


# GAMG example.  Requires OUTDomain,AIJ,1,<block_size>
feaprun-gamg:
	-@${MPIEXEC} -n $(NPROC) $(FEAPRUN) -ksp_type cg -ksp_monitor -log_summary -pc_type gamg -pc_gamg_type agg -pc_gamg_agg_nsmooths 1 -options_left

# Direct solvers.  These require OUTDomain,AIJ,<1>,<block_size> format
feaprun-superlu:
	-@${MPIEXEC} -n $(NPROC) $(FEAPRUN) -ksp_monitor_short -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package superlu_dist

feaprun-spooles:
	-@${MPIEXEC} -n $(NPROC) $(FEAPRUN) -ksp_monitor_short -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package spooles

feaprun-mumps:
	-@${MPIEXEC} -n $(NPROC) $(FEAPRUN) -ksp_monitor_short -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -memory_info -malloc_log -log_summary -on_error_attach_debugger

# Alternate form for job submission with differing protocal for defining number of processors and a machine file
# or
#	-@${MPIEXEC} -np $(NPROC) --machinefile machines $(FEAPRUN) -ksp_type cg -ksp_monitor -pc_type jacobi -log_summary -ksp_view -options_left

# Example submission line with bsub batch system and the prun command of the
# ETH-Zurich Brutus system  (WTIME is the exected maximum run time, and YESFILE
# is a file containing the inputs needed at the first prompt to accept the
# file names (i.e. a file with the single character 'y' followed by a carraige
# return).  Note that in batch submission your feapname file needs to be set-up already.
# this occurs automatically when you partition your problem with a single processor run.
feaprun-bsub:
	bsub -W $(WTIME) -n $(NPROC) -i $(YESFILE)  prun $(FEAPRUN) -ksp_type cg -ksp_monitor -pc_type jacobi -log_summary -ksp_view -options_left


# Note the make system sometimes complains about the following line but you can
# ignore the warning
.f.o:
	$(FF) -c $(FFOPTFLAG) -I$(FINCLUDE) $*.f -o $*.o

.c.o:
	$(CC) -c $(PETSC_CC_INCLUDES) $(CCOPTFLAG) $*.c -o $*.o

