c$Id:$
      subroutine pminvsqr(u,eq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'petscksp.h' include file                    24/07/2007
c       2. Remove obsolete petsc include files              07/03/2011
c       3. Only use processor's numpeq diagonal elments     07/01/2013
c          Remove numpeq argument
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute invers of square root of M-lumped

c      Inputs:
c        u(*)              - Temporary vector
c        eq(ndf,numnp)     - equation code array

c      Outputs:
c         Msqrt    - Square root of inverse mass
c         vneq     - Number of equations in eigenvectors (in "pfeapb.h")
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include  "finclude/petsc.h"
#     include  "pfeapb.h"
#     include  "pfeapc.h"

#     include  "cdata.h"
#     include  "fdata.h"
#     include  "iofile.h"
#     include  "sdata.h"
#     include  "setups.h"
#     include  "pointer.h"
#     include  "comblk.h"

      integer   i,j,k
      integer   eq(ndf,*)
      real*8    u(*)

c     Petsc values

      PetscScalar    sol_array(1)
      PetscOffset    i_sol
      PetscInt       starti,endi
      PetscErrorCode ierr
      PetscViewer    M_view

c     save

c     Determine vneq
c     These lines are not really correct the slightly over estimate the size
      vneq = 0
      do j = 1,numnp
        do i = 1,ndf
          if(eq(i,j).gt.0) then
            vneq = vneq + 1
          endif
        end do ! i
      end do ! j

c     Set global PETSc vectors

      if(Mdiag.eq.0) then
        write(iow,3000)
        call plstop()
      endif
      if(Msqrt.eq.0) then
        call VecDuplicate (Mdiag, Msqrt, ierr)
      endif

c     Set Mdiag values in u(1)

      call VecGetOwnershipRange(Mdiag, starti, endi, ierr)
      call VecGetArray         (Mdiag, sol_array, i_sol, ierr)

      do i = 1,numpeq
        u(i) = sol_array(i_sol+i)
        if(u(i).gt.0.0d0) then
          u(i) = 1.d0/sqrt(u(i))
        elseif(u(i).eq.0.0d0) then
          u(i) = 1.d+5
        elseif(u(i).lt.0.0d0) then
          write(iow,3001)
          call plstop()
        endif
      end do ! i

c     Insert values into Msqrt

      j = 0
      do i = 0,numpn*ndf-1
        k = mr(np(245)+i) - 1
        if(k.ge.0) then
          j = j + 1
          call VecSetValue     (Msqrt,k,u(j),INSERT_VALUES, ierr)
        endif
      end do ! i

      call VecAssemblyBegin(Msqrt, ierr)
      call VecAssemblyEnd  (Msqrt, ierr)

      call VecRestoreArray     (Mdiag, sol_array, i_sol, ierr)

      if(pfeap_view) then
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"msqr.m",M_view,
     &                            ierr)
        call PetscViewerSetFormat(M_view,PETSC_VIEWER_ASCII_MATLAB,
     &                            ierr)
        call VecView             (Msqrt,M_view,ierr)
        call PetscViewerDestroy  (M_view , ierr)
      endif


c     Formats

3000  format(' *ERROR* - Form lumped mass before solving eigenproblem')
3001  format(' *ERROR* - Negative lumped mass')
      end
