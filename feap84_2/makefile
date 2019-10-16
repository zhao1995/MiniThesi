
# N.B.  It is necessary to modify 'makefile.in' before using make.

include $(FEAPHOME8_4)/makefile.in

CLEANDIRS = include include/integer4 include/integer8 maintain \
	program contact elements plot unix user main windows window1 \
	window2 packages/arpack packages/arpack/archive packages/blas \
	packages/lapack packages/meshmod parfeap parfeap/packages/arpack \
	parfeap/partition program/memory unix/memory # parfeap/unix parfeap/windows 

feap: archive
	(cd main; make feap)
	@@echo "--> Feap executable made <--"
#	@@echo "--> N.B. Before parallel PETSc compile:"
#	@@echo "         Remove # from PETSc include files in:"
#	@@echo "         /parfeap/makefile and"
#	@@echo "         /parfeap/packages/arpack/makefile"

archive:   
	(cd program; make archive)
	(cd elements; make archive)
	(cd plot; make archive)
	(cd unix; make archive)
	(cd unix/nojpeg; make archive)   # Can be 'jpeg' instead
	(cd user; make archive)
	(cd unix/memory; make archive)
	(cd contact; make archive)
	@@echo "--> Feap Archive updated <--"

install: archive feap

clean:
	for i in $(CLEANDIRS); do (cd $$i; make clean); done
	if [ -f $(ARFEAP) ]; then rm $(ARFEAP); fi
	rcsclean -q
	@@echo "--> Feap rcscleaned <--"

#	(cd unix/memory; make archive)
#	(cd program/memory; make archive)
