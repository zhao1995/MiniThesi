c$Id:$
      subroutine paropk(n,u,v)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Author: D.S. Bindel 6/21/2005
c       Original version                                    01/11/2006
c       1. Add 'petscksp.h' include file                    24/07/2007
c       2. Remove obsolete petsc include files              07/03/2011
c       3. Change SetValue to loop over numpn*ndf eqs       07/03/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Apply K^-1 for Arnoldi or Lanczos

c      Inputs:
c         n    - Number of equations active (numpeq)
c         u(*) - Current vector

c      Outputs:
c         v(*) - Solution: K * v = u
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
#     include   "p_int.h"
#     include   "sdata.h"
#     include   "setups.h"
#     include   "pointer.h"
#     include   "comblk.h"

      integer    i,j,k, n
      real*8     u(*),v(*)

c     Petsc values

      PetscScalar    sol_array(1)
      PetscOffset    i_sol
      PetscInt       starti,endi
      PetscErrorCode ierr
      PetscViewer    Y_view
      PetscViewer    X_view
      PetscViewer    Z_view

c     save

c     Set global PETSc vectors

      if(mmode.eq.1) then

        if(xvec.eq.0) then
          call VecCreate        (PETSC_COMM_WORLD, xvec, ierr)
          call VecSetSizes      (xvec, numpeq, PETSC_DECIDE, ierr)
          call VecSetFromOptions(xvec, ierr)
          call VecDuplicate     (xvec, yvec, ierr)
          call VecDuplicate     (xvec, zvec, ierr)
        endif
        if(yvec.eq.0) then
          call VecDuplicate     (xvec, yvec, ierr)
        endif
        if(zvec.eq.0) then
          call VecDuplicate     (xvec, zvec, ierr)
        endif
        if(Msqrt.eq.0) then
          write(*,*) ' ERROR NO MASS SQRT'
          call plstop()
        endif

c       Multiply u(*) by Msqrt

        j = 0
        do i = 0,numpn*ndf-1
          k = mr(np(245)+i) - 1
          if(k.ge.0) then
            j = j + 1
            call VecSetValue   (xvec,k,u(j),INSERT_VALUES, ierr)
          endif
        end do ! i

        call VecAssemblyBegin(xvec, ierr)
        call VecAssemblyEnd  (xvec, ierr)

        call VecPointwiseMult  (yvec, Msqrt, xvec, ierr)

        if(pfeap_view) then
          call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"xvec.m",X_view,
     &                              ierr)
          call PetscViewerSetFormat(X_view,PETSC_VIEWER_ASCII_MATLAB,
     &                              ierr)
          call VecView             (yvec,X_view,ierr)
          call PetscViewerDestroy  (X_view , ierr)
        endif

        call MatMult           (Kmat, yvec, zvec, ierr)

        if(pfeap_view) then
          call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"zvec.m",Z_view,
     &                              ierr)
          call PetscViewerSetFormat(Z_view,PETSC_VIEWER_ASCII_MATLAB,
     &                              ierr)
          call VecView             (zvec,Z_view,ierr)
          call PetscViewerDestroy  (Z_view , ierr)
        endif

        call VecPointwiseMult  (yvec, Msqrt, zvec, ierr)

        if(pfeap_view) then
          call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"yvec.m",Y_view,
     &                              ierr)
          call PetscViewerSetFormat(Y_view,PETSC_VIEWER_ASCII_MATLAB,
     &                              ierr)
          call VecView             (yvec,Y_view,ierr)
          call PetscViewerDestroy  (Y_view , ierr)
        endif

        pfeap_view = .false.

c       Set my values in v(1)

        call VecGetOwnershipRange(yvec, starti, endi, ierr)
        call VecGetArray         (yvec, sol_array, i_sol, ierr)

        do i = 1,n
          v(i) = sol_array(i_sol+i)
        end do ! i

        call VecRestoreArray     (yvec, sol_array, i_sol, ierr)

      elseif(mmode.eq.3) then
        fp(1)  = na
        fp(2)  = nau
        fp(3)  = nal
        fp(4)  = np(20+npart)
        do i = 1,n
          v(i) = u(i)
        end do ! i
        call psolve(ittyp,v,fp,.false.,.true.,.true.,.false.)
      endif

      end
