!$Id:$
      subroutine umacr21(lct,ctl)

      USE OMP_LIB

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2015: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    01/11/2006
!       1. Remove 'prt' from argument list                  09/07/2009
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose:  User interface for adding solution command language
!                instructions.

!      Inputs:
!         lct       - Command character parameters
!         ctl(3)    - Command numerical parameters

!      Outputs:
!         N.B.  Users are responsible for command actions.  See
!               programmers manual for example.
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'compas.h'
      include  'iofile.h'
      include  'part0.h'
      include  'setups.h'
      include  'umac1.h'
      include  'pointer.h'
      include  'xpardiso.h'

      logical   pcomp,setvar,palloc,flags(5), pardfl
      character lct*15,tname*5
      real*8    ctl(3)

      save

!     Set command word

      if(pcomp(uct,'ma21',4)) then      ! Usual    form

        uct = 'pard'                    ! Specify 'pard'iso solver
        pardfl = .true.

      else                              ! Perform user operation

        npzs = -1                       ! Initialize number non-zeros


!       Delete any arrays if solver PARDISO set to OFF

        if(pcomp(lct,'off',3)) then

          if(np(npart).ne.0) then
            write(tname,'(4hTANG,i1)') npart
            setvar = palloc(npart,tname, 0,2)
          endif

          if(.not.solver) then
            flags(1:4) = .false.
            flags(  5) = .true.
            call usolve(flags,ctl)
          endif

          solver = .true.
          ittyp  =  -3               ! Default is profile

!         Output for deactivation

          if(ior.lt.0) then
            write(*,'(a)') ' --> PARDISO solver deactivated'
          endif

!       Solver properties for PARDISO

        else

!         Set 'solver' to false to indicate use of PARDISO solver

          solver = .false.
          ittyp  = -4                ! Prevent MFlop outputs

!         Set PARDISO control parameters

          nrhs        = 1            ! number right hand sides to solve
          maxfct      = 1            ! number with identical sparsity
          mnum        = 1            ! matrix number to use

!         Fill IPARM with values for solver

!-----[--.----+----.----+----.-----------------------------------------]
!         iparm( 1)    =  0           ! I: Filled with default values
          iparm( 1)    =  1           ! I: Furnishing values 2 to 64
          iparm( 2)    =  2           ! I: Fill-in reordering from METIS
!         iparm( 3)    =  nprocessor  ! I: Numbers of processors
          iparm( 3)    =  0           ! Reserved: set to zero
          iparm( 4)    =  0           ! I: No iterative-direct algorithm
          iparm( 5)    =  0           ! I: No fillin reduced permutation
          iparm( 6)    =  0           ! I: Solution x(1:neq)
          iparm( 7)    =  0           ! O: Number iterative refinements
          iparm( 8)    =  0           ! I: Maximum iterative refinements
          iparm( 9)    =  0           ! Reserved: set to zero
          iparm(10)    = 13           ! I: Perturb pivot by 1E-13
          iparm(11)    =  1           ! I: Nonsymmetric scaling
                                      !   (should be 0 for sym indef)
                                      !   (should be 1 for unsym)
          iparm(12)    =  0           ! I: Do not transpose A
          iparm(13)    =  0           ! I: Disable matching (for sym)
!         iparm(13)    =  1           ! I: Enable  matching (for unsym)
          iparm(14)    =  0           ! O: Number of perturbed pivots
          iparm(15)    =  0           ! O: Peak memory use
                                      !    (only in phase 1?)
          iparm(16)    =  0           ! O: Permanent memory use
                                      !    (only in phase 1?)
          iparm(17)    =  0           ! O: Size of factors and solution
          iparm(18)    =  1           ! I/O: Number nonzeros LDU (1=off)
          iparm(19)    =  1           ! I/O: Mflop factorization (1=off)
          iparm(20)    =  0           ! O: Numbers of CG Iterations
          iparm(21)    =  1           ! I: Pivoting 1x1 and 2x2
          iparm(22:64) =  0

!         Not previously set (IC = in-core; OOC = out of core)

!         iparm(60)    = 0            ! I: IC  mode
!         iparm(60)    = 1            ! I: IC  mode but OOC if big
!         iparm(60)    = 2            ! I: OOC mode

!         Initialize pt array to zero

          pt(1:64) = 0

!         Output of statistics

          if(nint(ctl(2)).eq.0) then
            istats = 0           ! No statistitcs output (DEFAULT)
          else
            istats = 1           ! Statistics output to screen
          endif 

!         Output for activation

          if(ior.lt.0 .or. pardfl) then
            write(*,'(a)') ' --> PARDISO solver activated'
            pardfl = .false.
          endif

        endif

      endif

      end
