c$Id:$
      subroutine usolve(flags,b)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2013: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Added MKL Pardiso to USSL                        12/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  User solver interface

c     Inputs:
c        flags(1) - Allocation and/or initialization phase
c        flags(2) - Perform factorization for direct solutions
c        flags(3) - Coefficient array unsymmetric
c        flags(4) - Solve equations
c        flags(5) - Purge any temporary storage used for solutions
c        b(*)     - RHS vector

c     Outputs:

c        flags(5) - True if error occurs (for solve/factor only)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include  'cdata.h'
      include  'compac.h'
      include  'debugs.h'
      include  'endata.h'
      include  'eqsym.h'
      include  'iofile.h'
      include  'part0.h'
      include  'setups.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'user_solvers.h'
      include 'mkl_pardiso.f77'

      character tname*5
      logical   flags(*), setvar,palloc
      integer   i,ldb,mrp, naux,nnu,ibsv
      real*8    b(*), dot, dparm(64)
      INTEGER idum(1)
      REAL*8  ddum(1)

      save

c     Set number right hand sides

      ldb  = neq
      mrp  = 0
      naux = 0
      nrhs = 1

c     Presolve setups

      if(flags(1)) then

       if(usol_pack.eq.1) then

!-------------------------------------
! MKL Pardiso Solver
!-------------------------------------

        maxfct = 1 !Maximum number of factors with identical nonzero sparsity structure
        mnum = 1 !Indicates the actual matrix for the solution phase
        mklerr = 0 ! initialize error flag

c       Check if problem is symmetric or unsymmetric

        if(flags(3)) then
          neqs = 1
          uall = .true.       ! Assemble total array
        else
          neqs = neq
          uall = .false.      ! Assemble upper triangular part
        endif

        ubycol = .false.      ! Store by rows
        udiag  = .true.       ! Diagonal in sparse array
        ulfl   = .false.      ! Prevents assembly of al in cassem

c       Establish current matrix indices

        call iters(nnu,-2)

c       Allocate matrix storage

        write(tname,'(4hTANG,i1)') npart
        setvar = palloc(npart,tname,nnu+1,2)
        call pzero(hr(np(npart)),nnu+1)

c       Adjust the pointer array

        setvar = palloc(227,'USOL3',neq+1,1)
        mr(np(227)) = 1
        do i = 1,neq
          mr(np(227)+i) = mr(np(225)+i-1) + 1
        end do ! i

            WRITE(iow,*) 'Number of nonzeros in factors = '
        if(uall) then ! Unsymmetric

          mtype = 1 ! real structurally symmetric
          !mtype = 11 ! real unsymmetric

C.. Initiliaze the internal solver memory pointer. This is only
C necessary for the FIRST call of the PARDISO solver.
        call PARDISOINIT(PT, MTYPE, IPARM)

C..
C.. Set up PARDISO control parameter
C..
          iparm(1) = 1 ! no solver default
          iparm(2) = 2 ! fill-in reordering from METIS
          iparm(3) = 1 ! numbers of processors
          iparm(4) = 0 ! no iterative-direct algorithm
          iparm(5) = 0 ! no user fill-in reducing permutation
          iparm(6) = 0 ! no copy solution into rhs
          iparm(7) = 0 ! not in use
          iparm(8) = 9 ! numbers of iterative refinement steps
          iparm(9) = 0 ! not in use
          iparm(10) = 13 ! perturbe the pivot elements with 1E-13
          iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
          iparm(12) = 0 ! not in use
          iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
          iparm(14) = 0 ! Output: number of perturbed pivots
          iparm(15) = 0 ! not in use
          iparm(16) = 0 ! not in use
          iparm(17) = 0 ! not in use
          iparm(18) = -1 ! Output: number of nonzeros in the factor LU
          iparm(19) = -1 ! Output: Mflops for LU factorization
          iparm(20) = 0 ! Output: Numbers of CG Iterations


        else ! symmetric

C.. Initiliaze the internal solver memory pointer. This is only
C necessary for the FIRST call of the PARDISO solver.
!      mtype = -2 by default, can be set to through variable 2 of FEAP macro command in umacr1.f
        call PARDISOINIT(PT, MTYPE, IPARM)

C..
C.. Set up PARDISO control parameter
C..
          iparm(1) = 1 ! no solver default
          iparm(2) = 2 ! fill-in reordering from METIS
          iparm(4) = 0 ! no iterative-direct algorithm
          iparm(5) = 0 ! no user fill-in reducing permutation
          iparm(6) = 0 ! no copy solution into rhs
          iparm(8) = 2 ! numbers of iterative refinement steps
          iparm(10) = 13 ! perturbe the pivot elements with 1E-13
          iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
          iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
          iparm(14) = 0 ! Output: number of perturbed pivots
          iparm(18) = -1 ! Output: number of nonzeros in the factor LU
          iparm(19) = -1 ! Output: Mflops for LU factorization
          iparm(20) = 0 ! Output: Numbers of CG Iterations
          mklerr = 0 ! initialize error flag

        endif !unsym/sym

       endif !usol_pack

c     Solution steps for assembled equations

      else

!-----------------------
c       Factor equations
!-----------------------

        if(flags(2)) then

         if(usol_pack.eq.1) then

!-------------------------------------
! MKL Pardiso Solver
!-------------------------------------

c         Ordering and Factorization

          !if(uall) then !unsymmetric

            phase = 12 ! reordering and factorization
            CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, 
     1           hr(np(npart)), mr(np(227)), mr(np(226)),
     2           idum, nrhs, iparm, msglvl, ddum, ddum, mklerr)

          !else !symmetric

         ! endif !unsym/sym


          WRITE(*,*) 'MKL - Reordering and factorization completed ... '
          IF (mklerr .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', mklerr
            call plstop()
          END IF


          if(debug) then
            WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
            WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)

            WRITE(iow,*) 'Number of nonzeros in factors = ',iparm(18)
            WRITE(iow,*) 'Number of factorization MFLOPS = ',iparm(19)
          endif

         endif !usol_pack

        endif !flag==2

!--------------------
c       Perform solve
!--------------------

        if(flags(4)) then

         if(usol_pack.eq.1) then

!-------------------------------------
! MKL Pardiso Solver
!-------------------------------------

          ! copy rhs for energy norm
          setvar = palloc(228,'USOL4',neq  ,2)
          ibsv   = np(228) - 1
          do i = 1,neq
            hr(ibsv+i) = b(i)
          end do ! i

          if(uall) then !unsymmetric

            iparm(8) = 2 ! max numbers of iterative refinement steps
            phase = 33 ! only factorization

            CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, 
     1      hr(np(npart)), mr(np(227)), mr(np(226)),
     2      idum, nrhs, iparm, msglvl, hr(np(228)), b, mklerr)

            WRITE(*,*) 'MKL - Solve completed ... '


          else !symmetric

            phase = 33 ! only factorization

            CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, 
     1      hr(np(npart)), mr(np(227)), mr(np(226)),
     2      idum, nrhs, iparm, msglvl, hr(np(228)), b, mklerr)

            WRITE(*,*) 'MKL - Solve completed ... '

          endif !unsym/sym


          if (mklerr .ne. 0) then
            write(*,*) 'The following ERROR was detected: ',mklerr
            call plstop()
          end if


c         Compute 'energy' for convergence test

          aengy  = dot(hr(np(228)),b,neq)
          setvar = palloc(228,'USOL4',    0,2)

         endif !usol_pack

        endif !flags==4

!--------------------------------
c       Purge storage in 'factor'
!--------------------------------

        if(flags(5)) then

         if(usol_pack.eq.1) then

!-------------------------------------
! MKL Pardiso Solver
!-------------------------------------

C.. Termination and release of memory
         phase = -1 ! release internal memory
         CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, ddum, idum,
     1    idum, idum, nrhs, iparm, msglvl, ddum, ddum, mklerr)

          setvar = palloc(227,'USOL3',    0,1)
          setvar = palloc(226,'USOL2',    0,1)  ! allocated in uiters
          setvar = palloc(225,'USOL1',    0,1)  ! allocated in uiters

         endif !usol_pack

        endif !flags==5

      endif !flags==1

      end
