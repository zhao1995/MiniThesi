c$Id:$
      subroutine FeapUserPCSetUp(Kmat,ksp,rhs,sol,ierr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'petscksp.h' include file                    24/07/2007
c       2. Change 'include/finclude' to 'finclude'          23/01/2009
c       3. Update petsc includes to v3.1                    20/07/2010
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Set up user defined PC interface for PETSc

c     The inputs to this routine are only illustrative.  They can and
c     probably should be changed in usolve.F depending upon the users
c     needs.  These routines will also need to store their own internal
c     data.

c     Inputs:
c        ksp      - ksp context
c        Kmat     - matrix
c        rhs      - right-hand side
c        sol      - solution

c     Outputs:

c        ierr     - PETSC error flag
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

#     include "finclude/petsc.h"

      Vec            rhs,sol
      Mat            Kmat
      KSP            ksp
      PetscErrorCode ierr

c     Common block to store data for user-provided preconditioner

      Vec              diag
      common /feapupc/ diag

#     include "iofile.h"

c     Example: Diagonal scaling stored in diag.  Could read in from
c              disk a matrix to store in common block /feapupc/

      call VecDuplicate(rhs,diag,ierr)
      call MatGetDiagonal(Kmat,diag,ierr)
      call VecReciprocal(diag,ierr)

c     Example: Initialize the solution.  Could read in from
c              disk an intial solution
      call VecSet(sol,1.d0,ierr)
      call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)

      write(iow,*) 'Setting up FEAP user PC'
      write(ilg,*) 'Setting up FEAP user PC'

      end

      subroutine FeapUserPCApply(dummy,x,y,ierr)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Apply user defined PC interface for PETSc

c     Inputs:
c        dummy    - optional user-defined contex, not used here
c        x        - input vector

c     Outputs:

c        y        - output vector = PC*x
c        ierr     - PETSC error flag (non-zero if there is an error)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

#     include "finclude/petsc.h"

      Vec            x,y
      PetscErrorCode ierr
      integer        dummy

c     Common block to store data for user-provided preconditioner

      Vec              diag
      common /feapupc/ diag

c     Apply diagonal preconditioner: y = diag .* x

      call VecPointwiseMult(y,x,diag,ierr)

      end
