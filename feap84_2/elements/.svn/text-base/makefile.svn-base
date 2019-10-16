include ../makefile.in

archive: 
	(cd couple3d; make archive)
	(cd frame; make archive)
	(cd material; make archive)
	(cd material/small; make archive)
	(cd material/finite; make archive)
	(cd shells; make archive)
	(cd solid1d; make archive)
	(cd solid2d; make archive)
	(cd solid3d; make archive)
	(cd thermal; make archive)
	(cd torsion; make archive)


clean: 
	(cd couple3d; make clean)
	(cd frame; make clean)
	(cd material; make clean)
	(cd material/small; make clean)
	(cd material/finite; make clean)
	(cd shells; make clean)
	(cd solid1d; make clean)
	(cd solid2d; make clean)
	(cd solid3d; make clean)
	(cd thermal; make clean)
	(cd torsion; make clean)
	rcsclean

checkout:
	co -q couple3d/makefile
	(cd couple3d; make checkout)
	co -q frame/makefile
	(cd frame; make checkout)
	co -q material/makefile
	(cd material; make checkout)
	co -q material/small/makefile
	(cd material/small; make checkout)
	co -q material/finite/makefile
	(cd material/finite; make checkout)
	co -q shells/makefile
	(cd shells; make checkout)
	co -q solid1d/makefile
	(cd solid1d; make checkout)
	co -q solid2d/makefile
	(cd solid2d; make checkout)
	co -q solid3d/makefile
	(cd solid3d; make checkout)
	co -q thermal/makefile
	(cd thermal; make checkout)
	co -q torsion/makefile
	(cd torsion; make checkout)

rcs:
	(cd couple3d; make rcs)
	(cd frame; make rcs)
	(cd material; make rcs)
	(cd material/small; make rcs)
	(cd material/finite; make rcs)
	(cd shells; make rcs)
	(cd solid1d; make rcs)
	(cd solid2d; make rcs)
	(cd solid3d; make rcs)
	(cd thermal; make rcs)
	(cd torsion; make rcs)
	ci -t-"" makefile

