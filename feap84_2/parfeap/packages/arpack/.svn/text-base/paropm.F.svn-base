c$Id:$
      subroutine paropm(n,u,v)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'petscksp.h' include file                    24/07/2007
c       2. Remove obsolete petsc include files              07/03/2011
c       3. Changed loop limts or sol_array read             05/01/2013
c          Change SetValue to loop over numpn*ndf eqs
c          Change input lengths to n
c       4. Count number of mass terms for subspace size     11/04/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute product of 'b' with 'v'

c      Inputs:
c         n        - Number of active equations (numpeq)
c         u(*)     - Iteration vectors

c      Outputs:
c         v(*)     - B * u
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include  "finclude/petsc.h"
#     include  "pfeapb.h"
#     include  "pfeapc.h"

#     include  "fdata.h"
#     include  "sdata.h"
#     include  "setups.h"
#     include  "pointer.h"
#     include  "comblk.h"

#     include  "parpacka.h"

      integer   i,j,k, n
      real*8    u(*),v(*)

      real*8    tma(1),tbuf(2)

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

c     Multiply u(*) by M

      j = 0
      do i = 0,numpn*ndf-1
        k = mr(np(245)+i) - 1
        if(k.ge.0) then
          j = j + 1
          call VecSetValue     (xvec,k,u(j),INSERT_VALUES, ierr)
        endif
      end do ! i

      call VecAssemblyBegin(xvec, ierr)
      call VecAssemblyEnd  (xvec, ierr)

      if(fl(1)) then
        call MatMult           (Mmat, xvec, yvec, ierr)
      else
        call VecPointwiseMult  (yvec, Mdiag, xvec, ierr)
      endif

c     Get my values in v(i)

      call VecGetOwnershipRange(yvec, starti, endi, ierr)
      call VecGetArray         (yvec, sol_array, i_sol, ierr)

      do i = 1,n
        v(i) = sol_array(i_sol+i)
      end do ! i

c     Count number of non-zero mass to set subspace size

      if(parmassfl) then
        nummass = 0
        do i = 1,n
          if(v(i).gt.0.0d0) nummass = nummass + 1
        end do ! i
        tma(1) = dble(nummass)
        call pfeapsr(tma,tbuf,1,.true.)
        nummasst = nint(tma(1))
      else
        call VecRestoreArray     (yvec, sol_array, i_sol, ierr)
      endif

      end
