c$Id:$
      subroutine ppmodin(phi, y, u,v,tu,nv,tneq,mtyp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Added petscksp.h                                 01/05/2007
c       2. Change 'include/finclude' to 'finclude'          23/01/2009
c       3. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c       4. Update petsc includes to v3.1                    20/07/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Initial Conditions for modal solution of
c               linear transient problems

c      Inputs:
c         phi(vneq,*) - Mass orthonormal eigenvectors
c         u(*)        - Initial displacements
c         v(*)        - Initial velocities
c         nv          - Number of eigenpairs (converged)
c         tneq        - Number of d.o.f. in model
c         mtyp        - Type of mass (true for consistent,false for lumped)

c      Scratch:
c         tu(*,2)     - Temporary vector

c      Outputs:
c         y(nv,3)     - Eigensolution at time t_0
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include  "finclude/petsc.h"
#     include  "pfeapc.h"

#     include  "pfeapb.h"

#     include  "sdata.h"
#     include  "setups.h"
#     include  "pointer.h"
#     include  "comblk.h"

      logical   mtyp
      integer   i,j,k,l,n, nv,tneq
      real*8    phi(vneq,*), y(nv,3), u(*),v(*), tu(tneq,*)

c     Petsc values

      PetscBool      chk
      PetscErrorCode ierr

      save

c     Compute 'z' and 'b' projection to form 'h'

      if(xvec.eq.0) then
        call VecCreate        (PETSC_COMM_WORLD, xvec, ierr)
        call VecSetSizes      (xvec, numpeq, PETSC_DECIDE, ierr)
        call VecSetFromOptions(xvec, ierr)
      endif

      if(yvec.eq.0) then
        call VecDuplicate     (xvec, yvec, ierr)
      endif

c     Multiply v(*,j) by M

      do j = 1,2
        l = 0
        do i = 0,numpn*ndf-1
          k = mr(np(245)+i) - 1
          if(k.ge.0) then
            l = l + 1
            call VecSetValue(xvec,k,tu(l,j),INSERT_VALUES, ierr)
          endif
        end do ! i

        call VecAssemblyBegin(xvec, ierr)
        call VecAssemblyEnd  (xvec, ierr)

        if(mtyp) then
          call MatMult( Mmat, xvec, yvec, ierr )
        else
          call VecPointwiseMult( yvec, Mdiag, xvec, ierr )
        endif

        do n = 1,nv
          l = 0
          do i = 0,numpn*ndf-1
            k = mr(np(245)+i) - 1
            if(k.ge.0) then
              l = l + 1
              call VecSetValue(xvec,k,phi(l,n),INSERT_VALUES, ierr)
            endif
          end do ! i

          call VecAssemblyBegin(xvec, ierr)
          call VecAssemblyEnd  (xvec, ierr)

          call VecDot(xvec,yvec, y(n,j), ierr)
        end do ! n
      end do ! j

      end
