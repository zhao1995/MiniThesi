# fortran compiler = intel fortran
fcomp = ifort
switch = -module obj -w -f77rtl -fpp -D_BATCH_ -D_NOLIC_ -D_LINUX_ -D__INTEL_ -D__SHMDEBUG_ -O3 -qopenmp -parallel -intconstant -zero -assume:byterecl -check:pointers -check:bounds -check:uninit
#switch = -module obj2 -w -f77rtl -fpp -D_BATCH_ -D_NOLIC_ -D_LINUX_ -D__INTEL_ -D__SHMDEBUG_ -O3 -intconstant -zero -assume:byterecl -check:pointers -check:bounds -check:uninit

baseList = src/arclen.for src/bfgs.for src/contact.for src/cor.for src/crit.for src/csr_help.for \
src/curve.for src/cvuser.for src/cylt1.for src/cylt2.for src/cylt3.for src/cyltalg.for \
src/dynamics.for src/eig1.for src/ext.for src/feap.for src/feap_micro.for \
src/feaps1.for src/feaps2.for src/feast.for src/fepost.for \
src/help.for src/hpgl.for src/iniuser.for src/lanczos.for src/matelib3d.for \
src/mat.for src/pardisoint.for src/pbcg.for src/pgmres2.for src/pgmres.for \
src/pmesh1.for src/precond.for \
src/reme_nl1.for src/reme_nl2.for src/remesh1.for src/remesh2.for src/remesh3.for  \
src/simplex_qld.for src/smoo_control.for src/smsolv.for src/subsp.for src/superluint.for  \
src/ueigen.for \
src/unika.for src/section.for src/dplot.for src/tecplot.for src/paraview.for \
src_int/plot_intel.for src_int/plot_intel_cw.for src_int/feap_shmem.for
baseObj = $(baseList:.for=.o)

eleList = src/elmt01s.for src/elmt02s.for src/elmt03s.for src/elmt04s.for src/elmt05s.for src/elmt06s.for \
src/elmt07s.for src/elmt08s.for src/elmt09s.for src/elmt10s.for src/elmt11s.for src/elmt12s.for \
src/elmt13s.for src/elmt14s.for src/elmt15s.for src/elmt16s.for src/elmt17s.for src/elmt18s.for \
src/elmt19s.for src/elmt20s.for src/elmt21s.for src/elmt22s.for src/elmt23s.for src/elmt24s.for \
src/elmt25s.for src/elmt26s.for src/elmt29s.for src/elmt30s.for src/elmt31s.for src/elmt32s.for \
src/elmt33s.for src/elmt34s.for src/elmt45s_sshell_eas.for src_rwth/elmt61s.for src_rwth/elmt62s.for \
src_rwth/elmt63s.for
eleObj = $(eleList:.for=.o)

matList = src/mate3d01.for src/mate3d02.for src/mate3d03.for src/mate3d04.for src/mate3d05.for \
src/mate3d06.for src/mate3d07.for src/mate3d08.for src/mate3d09.for src/mate3d10.for src/mate3d11.for \
src/mate3d12.for src/mate3d14.for src/mate3d15.for src_rwth/mate3d16.for src_rwth/mate3d17.for
matObj = $(matList:.for=.o)

modList = src_rwth/module.for src_int/allocation_int.for
modObj = $(modList:.for=.o)

rwthList = src_rwth/dalloc.for src_rwth/feaps3.for src_rwth/feaps4.for src_rwth/feaps5.for src_rwth/feaps6.for \
src_rwth/elmlib.for src_rwth/isogeo.for  \
src_rwth/rhinoplot.for src_rwth/plot.for \
src_rwth/pplotf.for src_rwth/gframe.for src_rwth/pmacr.for src_rwth/pofeap.for
rwthObj = $(rwthList:.for=.o)

all: modules base ele mat feap

base: $(baseObj)
%.o: %.for
	$(fcomp) $(switch) -c $< -o obj/$(notdir $@)

ele: $(eleObj)
%.o: %.for
	$(fcomp) $(switch) -c $< -o obj/$(notdir $@)

mat: $(matObj)
%.o: %.for
	$(fcomp) $(switch) -c $< -o obj/$(notdir $@)

modules: $(modObj)
%.o: %.for
	$(fcomp) $(switch) -c $< -o obj/$(notdir $@)

rwth: $(rwthObj)
%.o: %.for
	$(fcomp) $(switch) -c $< -o obj/$(notdir $@)

feap: modules base ele mat rwth
	$(fcomp) $(switch) obj/*.o -o obj/feap_console -mkl -lrt
clean:
	@echo "target clean: delete obj/*"
	rm -f obj/*
