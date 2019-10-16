include ../makefile.in

archive:
	(cd main;  make archive)
	(cd ntrnd; make archive)
	(cd nts2d; make archive)
	(cd nts3d; make archive)
	(cd ptpnd; make archive)
	(cd tie2d; make archive)
	(cd util;  make archive)

clean:
	(cd main; make clean)
	(cd ntrnd; make clean)
	(cd nts2d; make clean)
	(cd nts3d; make clean)
	(cd ptpnd; make clean)
	(cd tie2d; make clean)
	(cd util; make clean)
	rcsclean

checkout:
	co main/makefile
	(cd main; make checkout)
	co ntrnd/makefile
	(cd ntrnd; make checkout)
	co nts2d/makefile
	(cd nts2d; make checkout)
	co nts3d/makefile
	(cd nts3d; make checkout)
	co ptpnd/makefile
	(cd ptpnd; make checkout)
	co tie2d/makefile
	(cd tie2d; make checkout)
	co util/makefile
	(cd util; make checkout)

rcs:
	(cd main; make rcs)
	(cd ntrnd; make rcs)
	(cd nts2d; make rcs)
	(cd nts3d; make rcs)
	(cd ptpnd; make rcs)
	(cd tie2d; make rcs)
	(cd util; make rcs)
	ci -t-"" makefile

