c$Id:$
      subroutine upremas(fl)

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
c       5. Change 'ndf' to 'nbsk' to match usolve.F         05/01/2013
c       6. Update for loss of VecValid, MatValid            05/01/2013
c       7. Update matrix creation calls                     06/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Mass interface for PETSc

c     Inputs:
c        fl(1) - Form Consistent mass if true
c        fl(2) - Form Lumped mass if true

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "finclude/petsc.h"
#     include   "pfeapc.h"

#     include   "cdata.h"
#     include   "comblk.h"
#     include   "endata.h"
#     include   "iofile.h"
#     include   "pointer.h"
#     include   "sdata.h"
#     include   "setups.h"
#     include   "pfeapb.h"

      PetscErrorCode ierr
      PetscBool      chk

      logical    fl(2)

c     save

c     Preform setup

      if(fl(1)) then   ! Consistent mass allocate
        if(Mmat.eq.0) then
          call MatCreate(PETSC_COMM_WORLD,Mmat,ierr)
          call MatSetSizes(Mmat,numpeq,numpeq,PETSC_DETERMINE,
     &                        PETSC_DETERMINE,ierr)
          if(pfeap_bcin) call MatSetBlockSize(Mmat,nsbk,ierr)
          call MatSetFromOptions(Mmat, ierr)
          if(pfeap_blk) then
            call MatSetType(Mmat,MATBAIJ,ierr)
            call MatMPIBAIJSetPreallocation(Mmat,nsbk,
     &                               PETSC_NULL_INTEGER,mr(np(246)),
     &                               PETSC_NULL_INTEGER,mr(np(247)),
     &                               ierr)
            call MatSeqBAIJSetPreallocation(Mmat,nsbk,
     &                               PETSC_NULL_INTEGER,mr(np(246)),
     &                               ierr)
          else
            call MatSetType(Mmat,MATAIJ,ierr)
            call MatMPIAIJSetPreallocation(Mmat,
     &                               PETSC_NULL_INTEGER,mr(np(246)),
     &                               PETSC_NULL_INTEGER,mr(np(247)),
     &                               ierr)
            call MatSeqAIJSetPreallocation(Mmat,
     &                               PETSC_NULL_INTEGER,mr(np(246)),
     &                               ierr)
          endif
        endif
      elseif(fl(2)) then   ! Lumped mass allocate

        if(Mdiag.eq.0) then
          call VecCreate        (PETSC_COMM_WORLD, Mdiag, ierr)
          call VecSetSizes      (Mdiag, numpeq, PETSC_DECIDE, ierr)
          call VecSetFromOptions(Mdiag, ierr)
        endif
      elseif(.not.fl(1) .and. .not.fl(2)) then
        write(iow,*) ' ERROR DID NOT ALLOCATE MASS MATRIX'
        write(  *,*) ' ERROR DID NOT ALLOCATE MASS MATRIX'
        call plstop()
      endif
      if(ierr .ne. 0) then
        write(iow,*) 'Error on MatCreate'
        write(  *,*) 'Error on MatCreate'
        call plstop()
      endif

      if(fl(1)) then
        call MatZeroEntries   (Mmat,ierr)
      else
        call VecZeroEntries   (Mdiag,ierr)
      endif

c     Vector for matrix multiply

      if(xvec.eq.0) then
        call VecCreate        (PETSC_COMM_WORLD, xvec, ierr)
        call VecSetSizes      (xvec, numpeq, PETSC_DECIDE, ierr)
        call VecSetFromOptions(xvec, ierr)
      endif

      end
