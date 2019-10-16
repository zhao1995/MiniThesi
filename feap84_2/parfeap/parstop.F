c$Id:$
      subroutine parstop()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add a ksp object destroy                         01/05/2007
c       2. Add 'petscksp.h' include file                    24/07/2007
c       3. Change 'include/finclude' to 'finclude'          23/01/2009
c       4. Update petsc includes to v3.1                    20/07/2010
c       5. Update for loss of VecValid, MatValid            05/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Close any open parallel array and delete memory use
c               PETSc version

c      Inputs:
c         none

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "finclude/petsc.h"
#     include   "pfeapb.h"
#     include   "pfeapc.h"

      PetscErrorCode ierr
      PetscBool      chk

      save

c     Close out PETSc matrices

      if(Kmat.ne.0) then
        call MatDestroy(Kmat, ierr)
      endif

      if(Pmat.ne.0) then
        call MatDestroy(Pmat, ierr)
      endif

      if(Mmat.ne.0) then
        call MatDestroy(Mmat, ierr)
      endif

      if(pfeap_dstr) then
        call KSPDestroy(kspsol, ierr)
      endif

c     Close out PETSc vectors


      if(Mdiag.ne.0) then
        call VecDestroy(Mdiag, ierr)
      endif

      if(Msqrt.ne.0) then
        call VecDestroy(Msqrt, ierr)
      endif

      if(xvec.ne.0) then
        call VecDestroy(xvec, ierr)
      endif

      if(yvec.ne.0) then
        call VecDestroy(yvec, ierr)
      endif

      if(zvec.ne.0) then
        call VecDestroy(zvec, ierr)
      endif

      if(rhs.ne.0) then
        call VecDestroy(rhs , ierr)
      endif

      if(sol.ne.0) then
        call VecDestroy(sol , ierr)
      endif

c     Stop PETSc

      call PetscFinalize(ierr)

      end
