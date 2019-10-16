c$Id:$
      subroutine parbmat(phib,mb,cmass,lmass,t, neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Added petscksp.h                                 01/05/2007
c       2. Change 'include/finclude' to 'finclude'          23/01/2009
c       3. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c       4. Update petsc includes to v3.1                    07/20/2010
c       5. Update for removal of VecValid                   05/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute base matrix coupling term

c     Inputs:

c     Outputs:

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include   "finclude/petsc.h"

#     include   "fdata.h"
#     include   "pfeapb.h"
#     include   "pfeapc.h"
#     include   "sdata.h"
#     include   "setups.h"
#     include   "pointer.h"
#     include   "comblk.h"

c     Petsc values

      PetscScalar    sol_array(1)
      PetscOffset    i_sol
      PetscInt       starti,endi
      PetscBool      chk
      PetscErrorCode ierr

      integer    neq, i,j,k
      real*8     t(*),cmass(*),lmass(*),phib(*),mb(*)

      save

c     Parallel solution

c     Set PETSc vectors

      if(xvec .eq. 0) then
        call VecCreate        (PETSC_COMM_WORLD, xvec, ierr)
        call VecSetSizes      (xvec, numpeq, PETSC_DECIDE, ierr)
        call VecSetFromOptions(xvec, ierr)
      endif

      if(yvec .eq. 0) then
        call VecDuplicate     (xvec, yvec, ierr)
      endif

c     Multiply phib(*) by M

      j = 0
      do i = 0,ndf*numpn-1
        k = mr(np(245)+i) - 1
        if(k.ge.0) then
          j = j + 1
          call VecSetValue(xvec,k,phib(j),INSERT_VALUES, ierr)
        endif
      end do ! i
      call VecAssemblyBegin(xvec, ierr)
      call VecAssemblyEnd  (xvec, ierr)

      if(fl(1)) then
        call MatMult( Mmat, xvec, yvec, ierr )
      else
        call VecPointwiseMult( yvec, Mdiag, xvec, ierr )
      endif

c     Set my values in t(*)

      call VecGetOwnershipRange(yvec, starti, endi, ierr)
      call VecGetArray         (yvec, sol_array, i_sol, ierr)

      do i = 1,numpeq
        t(i) = sol_array(i_sol+i) - mb(i)
      end do ! i

      call VecRestoreArray     (yvec, sol_array, i_sol, ierr)

      end
