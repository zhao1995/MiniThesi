c$Id:$
      subroutine parmv(n,v)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Author: D.S. Bindel 6/21/2005
c       Original version                                    01/11/2006
c       1. Add 'petscksp.h' include file                    24/07/2007
c       2. Remove obsolete petsc include files              07/03/2011
c       3. Adjust SetValue to loop over numpn*ndf eqs       07/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute eigenvectors for mode 1 arpack solution using
c               lumped mass.

c      Inputs:
c         n    - Number of equations active (numpeq)
c         v(*) - Current Ritz vector

c      Outputs:
c         v(*) - Matrix multiply: Msqrt * v -> v
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include  "finclude/petsc.h"
#     include  "pfeapb.h"
#     include  "pfeapc.h"

#     include   "cdata.h"
#     include   "compas.h"
#     include   "evdata.h"
#     include   "fdata.h"
#     include   "ndata.h"
#     include   "part0.h"
#     include   "sdata.h"
#     include   "setups.h"
#     include   "pointer.h"
#     include   "comblk.h"

      integer    i,j,k, n
      real*8     v(*)

c     Petsc values

      PetscScalar    sol_array(1)
      PetscOffset    i_sol
      PetscInt       starti,endi
      PetscErrorCode ierr

c     save

c     Set global PETSc vectors

      if(xvec.eq.0) then
        call VecCreate        (PETSC_COMM_WORLD, xvec, ierr)
        call VecSetSizes      (xvec, numpeq, PETSC_DECIDE, ierr)
        call VecSetFromOptions(xvec, ierr)
        call VecDuplicate     (xvec, yvec, ierr)
      endif
      if(yvec.eq.0) then
        call VecDuplicate     (xvec, yvec, ierr)
      endif

c     Multiply v(*) by Msqt

      j = 0
      do i = 0,numpn*ndf-1
        k = mr(np(245)+i) - 1
        if(k.ge.0) then
          j = j + 1
          call VecSetValue   (xvec,k,v(j),INSERT_VALUES, ierr)
        endif
      end do ! i

      call VecAssemblyBegin(xvec, ierr)
      call VecAssemblyEnd  (xvec, ierr)

      call VecPointwiseMult  (yvec, Msqrt, xvec, ierr)

c     Set my values in v(1)

      call VecGetOwnershipRange(yvec, starti, endi, ierr)
      call VecGetArray         (yvec, sol_array, i_sol, ierr)

      do i = 1,n
        v(i) = sol_array(i_sol+i)
      end do ! i

      call VecRestoreArray     (yvec, sol_array, i_sol, ierr)

      end
