
! Common variables

      integer usol_pack

      common /user_solvers/ usol_pack

! MKL Pardiso variables

      INTEGER*8 pt(64)
      INTEGER maxfct, mnum, mtype, phase, nrhs, mklerr, msglvl
      INTEGER iparm(64)

      common /user_solvers/ maxfct, mnum, mtype, 
     &  phase, nrhs, mklerr, msglvl, iparm


      common /user_solvers/ pt
