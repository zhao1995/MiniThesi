c$Id:$
      subroutine usolve(flags,b)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise send/receive (pmaxr,pmaxs use)            24/11/2006
c       2. Correct computation of pmaxr and pmaxs           16/03/2007
c       3. Move kspsol to pfeapc.h                          01/05/2007
c       4. Move KSPSetFromOptions before set precond.       24/07/2007
c          Add 'petscksp.h' include file
c       5. Change 'include/finclude' to 'finclude'          23/01/2009
c       6. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c       7. Add flush after outputs to force display of text 06/02/2009
c       8. Add call to pflush                               25/02/2009
c       9. Add allocation of at least one TEMP0 word for    24/11/2009
c          one processor runs
c      10. Update petsc includes to v3.1                    20/07/2010
c      11. Update petsc call to v3.3                        20/12/2012
c      12. Update for loss of VecValid, MatValid            05/01/2013
c      13. Update matrix allocation calls                   06/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Solver interface for PETSc

c     Inputs:
c        flags(1) - Allocation and/or initialization phase
c        flags(2) - Perform factorization for direct solutions
c        flags(3) - Coefficient array unsymmetric
c        flags(4) - Solve equations
c        flags(5) - Purge storage of pointers
c        b(*)     - RHS vector

c     Outputs:

c        flags(5) - True if error occurs (for factor/solve only)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "finclude/petsc.h"
#     include   "pfeapb.h"
#     include   "pfeapc.h"

#     include   "cdata.h"
#     include   "comblk.h"
#     include   "compas.h"
#     include   "debugs.h"
#     include   "endata.h"
#     include   "idptr.h"
#     include   "iofile.h"
#     include   "pointer.h"
#     include   "print.h"
#     include   "sdata.h"
#     include   "setups.h"

      PetscScalar    sol_array(1)
      PetscOffset    i_sol
      PetscInt       starti,endi, iters
      PetscErrorCode ierr
      PetscBool      chk
      PetscViewer    K_view, B_view
      KSPConvergedReason reason
      PC             pc
      PetscReal      p_rtol,p_atol,p_dtol
      PetscInt       p_maxit

      real*8         info(MAT_INFO_SIZE), nza,nzr,mal

      logical    flags(5), setvar, palloc
      character  creason(7)*29,nreason(10)*26, pctype*20

      external   FeapUserPCApply

      integer    i, j, k, pmaxr, pmaxs
      real*8     b(*)

c     save

      data       creason / ' Unknown -1',
     &                     ' Satisfied residual tolerance',
     &                     ' Satisfied absolute tolerance',
     &                     ' Satisfied iteration limit',
     &                     ' Unused - 5',
     &                     ' Unused - 6',
     &                     ' Satisfied step length'/

      data       nreason / ' Unknown ',
     &                     ' Solution null',
     &                     ' Iterations exceeded',
     &                     ' Solution diverged',
     &                     ' Breakdown?',
     &                     ' Breakdown BiCG?',
     &                     ' Matrix nonsymmetric',
     &                     ' Indefinite Preconditioner',
     &                     ' NaN in solution',
     &                     ' Indefinite Matrix'/

c     Presolve setups

      if(flags(1)) then

        if(Kmat.eq.0) then
          call MatCreate(PETSC_COMM_WORLD,Kmat,ierr)
          call MatSetSizes(Kmat,numpeq,numpeq,PETSC_DETERMINE,
     &                        PETSC_DETERMINE,ierr)
          if(pfeap_bcin) call MatSetBlockSize(Kmat,nsbk,ierr)
          call MatSetFromOptions(Kmat, ierr)
          if(pfeap_blk) then
            call MatSetType(Kmat,MATBAIJ,ierr)
            call MatMPIBAIJSetPreallocation(Kmat,nsbk,
     &                               PETSC_NULL_INTEGER,mr(np(246)),
     &                               PETSC_NULL_INTEGER,mr(np(247)),
     &                               ierr)
            call MatSeqBAIJSetPreallocation(Kmat,nsbk,
     &                               PETSC_NULL_INTEGER,mr(np(246)),
     &                               ierr)
          else
            call MatSetType(Kmat,MATAIJ,ierr)
            call MatMPIAIJSetPreallocation(Kmat,
     &                               PETSC_NULL_INTEGER,mr(np(246)),
     &                               PETSC_NULL_INTEGER,mr(np(247)),
     &                               ierr)
            call MatSeqAIJSetPreallocation(Kmat,
     &                               PETSC_NULL_INTEGER,mr(np(246)),
     &                               ierr)
          endif
          if(ierr .ne. 0) then
            flags(5) = .true.
            write(*,*) 'Error on MatCreate'
          endif
        endif
        call MatZeroEntries     (Kmat,ierr)

c     Solution steps for assembled equations

      else

c       Factor equations

        if(flags(2)) then

c         Complete assembly of stiffness matrix

          call MatAssemblyBegin (Kmat, MAT_FINAL_ASSEMBLY, ierr)
          call MatAssemblyEnd   (Kmat, MAT_FINAL_ASSEMBLY, ierr)
          call MatSetOption(Kmat,MAT_NEW_NONZERO_LOCATIONS,
     &                      PETSC_TRUE,ierr)

c         Check on preallocation of matrix

          if(debug) then
            call MatGetInfo(Kmat, MAT_LOCAL, info, ierr)
            nza = info(MAT_INFO_NZ_ALLOCATED)
            nzr = info(MAT_INFO_NZ_USED)
            mal = info(MAT_INFO_MALLOCS)
            write(iow,2000) ' Matrix Information',rank,nza,nzr,mal
            write(  *,2000) ' Matrix Information',rank,nza,nzr,mal
2000  format(/a/' Rank[',i3,']'/
     &          ' Matrix non-zeros allocated =',1p,1d21.12/
     &          ' Matrix non-zeros used      =',1p,1d21.12/
     &          ' Number of mallocs          =',1p,1d21.12)

          endif

c         Stiffness matrix viewer

          if(pfeap_view) then
            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"stiff.m",K_view,
     &                                ierr)
            call PetscViewerSetFormat(K_view,PETSC_VIEWER_ASCII_MATLAB,
     &                                ierr)
            call MatView             (Kmat,K_view,ierr)
            call PetscViewerDestroy  (K_view,ierr)
          endif

        endif

c       Perform assembly of RHS

        if(flags(4)) then

          if(rhs.eq.0) then
            call VecCreate        (PETSC_COMM_WORLD, rhs, ierr)
            call VecSetSizes      (rhs, numpeq, PETSC_DECIDE, ierr)
            call VecSetFromOptions(rhs, ierr)
          endif

          if(sol.eq.0) then
            call VecDuplicate     (rhs, sol, ierr)
          endif

          j = 0
          do i = 0, ndf*numpn-1
            k = mr(np(245) + i)
            if( k .gt. 0 ) then
              j = j + 1
              call VecSetValue  (rhs, k-1, b(j), INSERT_VALUES, ierr)
            endif
          end do

          call VecAssemblyBegin (rhs, ierr)
          call VecAssemblyEnd   (rhs, ierr)

c DEBUGGING
c         call VecGetOwnershipRange(rhs, starti, endi ,ierr)
c         call VecGetArray         (rhs, sol_array, i_sol, ierr)
c         call mprint(sol_array(i_sol),1,numpeq,1,'RHS-USOLVE')
c         call VecRestoreArray(rhs, sol_array, i_sol, ierr)
c DEBUGGING


c         RHS Vector viewer

          if(pfeap_view) then
            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"rhs.m",B_view,
     &                                ierr)
            call PetscViewerSetFormat(B_view,PETSC_VIEWER_ASCII_MATLAB,
     &                                ierr)
            call VecView             (rhs,B_view,ierr)
            call PetscViewerDestroy  (B_view,ierr)
          endif

c         Create solution handle solution

          if(pfeap_tang) then
            if(pfeap_dstr) then
              call KSPDestroy       (kspsol,           ierr)
            endif
            pfeap_dstr = .true.
            call KSPCreate        (PETSC_COMM_WORLD, kspsol   ,ierr)
            call KSPSetOperators  (kspsol, Kmat, Kmat,
     &                             DIFFERENT_NONZERO_PATTERN  ,ierr)

c           Options from command line

            call KSPSetFromOptions(kspsol,ierr)

c           Set preconditioner type

            call KSPGetPC             (kspsol, pc ,      ierr)
            call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-pc_type',
     &           pctype,chk,ierr)

            if (chk .eqv. PETSC_TRUE) then

              if(pctype.eq.'shell') then

                call PCShellSetApply(pc,FeapUserPCApply,ierr)
                call FeapUserPCSetUp(Kmat,kspsol,rhs,sol,ierr)

              else if   (pctype.eq.'gamg'
     &               .or.pctype.eq.'prometheus'
     &               .or.pctype.eq.'ml') then

                call PCSetCoordinates(pc,ndm,numpn,hr(np(43)),ierr)

              endif

            endif
            pfeap_tang = .false.
          endif

c         Set tolerance values for convergence checks

          call KSPGetTolerances (kspsol,p_rtol,p_atol,p_dtol,p_maxit,
     &                           ierr)
          if(pmaxit.le.0) then
            pmaxit = p_maxit
          endif
          call KSPSetTolerances (kspsol,  itol,  atol,  dtol,pmaxit,
     &                           ierr)                   ! From command data

c         Solve equations and destroy handle

          call KSPSolve         (kspsol, rhs, sol, ierr)
          call KSPGetConvergedReason(kspsol,reason,ierr)
          if(reason.gt.0 .and. echo) then
            call KSPGetIterationNumber(kspsol,iters,ierr)
            write(  *,*) ' CONVERGENCE: ',creason(reason),
     &                   ' Iterations = ',iters
            write(iow,*) ' CONVERGENCE: ',creason(reason),
     &                   ' Iterations = ',iters
          elseif(reason.lt.0) then
            write(  *,*) ' NO CONVERGENCE REASON: ',nreason(-reason)
            write(iow,*) ' NO CONVERGENCE REASON: ',nreason(-reason)
          endif
c         call KSPDestroy       (kspsol,           ierr)

c         Set local values in b(*)

          call VecGetOwnershipRange(sol, starti, endi ,ierr)
          call VecGetArray         (sol, sol_array, i_sol, ierr)

          do i = 1,numpeq
            b(i) = sol_array(i_sol + i)
          end do ! i
c         call mprint(b,1,numpeq,1,'SOL-USOLVE')

          call VecRestoreArray(sol, sol_array, i_sol, ierr)

c         Compute max message buffer sizes for ghost communication

          pmaxr = 0
          pmaxs = 0
          do i = 0, ntasks-1
            pmaxr = max( pmaxr, mr(np(248)+i) )
            pmaxs = pmaxs + mr(np(250)+i)
          end do
          pmaxr = pmaxr*ndf
          pmaxs = pmaxs*ndf

c         Perform communication and set of ghost values in solution

          setvar = palloc( 120, 'TEMP0', max(1,pmaxr+pmaxs), 2)
          call psetb(b,mr(np(248)),mr(np(249)),mr(np(250)),mr(np(251)),
     &                 mr(id31),ndf, hr(np(120)),hr(np(120)+pmaxr))
          setvar = palloc( 120, 'TEMP0',    0, 2)


c         Compute 'energy' for convergence test

          call VecDot(rhs,sol, aengy,ierr)

        endif

c       Purge storage in 'factor'

        if(flags(5)) then
c         write(*,*) 'Rank: ',rank,' USOLVE flag(5) true'
        endif

      endif

c     if(rank.eq.0) call pflush(-1)
      call pflush(-1)

      end
