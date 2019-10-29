!$Id:$
      subroutine usolve(flags,b)

      use OMP_LIB

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2015: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    01/11/2006
!-----[--.----+----.----+----.-----------------------------------------]
!     Purpose:  Interface to Intel Pardiso MKL solver

!     Inputs:
!        flags(1) - Allocation and/or initialization phase
!        flags(2) - Perform factorization for direct solutions
!        flags(3) - Coefficient array unsymmetric
!        flags(4) - Solve equations
!        flags(5) - Purge any temporary storage used for solutions
!        b(*)     - RHS vector

!     Outputs:

!        flags(5) - True if error occurs (for solve/factor only)
!-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'compac.h'
      include   'debugs.h'
      include   'endata.h'
      include   'iofile.h'
      include   'part0.h'
      include   'setups.h'

      include   'pointer.h'
      include   'comblk.h'

      include   'xpardiso.h'

      character        :: tname*5
      logical          :: flags(*), setvar,palloc
      real    (kind=8) :: b(*),ddum
!     integer (kind=8) :: pt(64)
      integer (kind=4) :: phase, error, idum
      integer (kind=4) :: i, npzsv

      save

!     Presolve setups

      if(flags(1)) then

!       Nonsymmetric solver

        if(flags(3)) then

          ubycol = .false.    ! Sparse storage by columns
          udiag  = .true.     ! Store diagonal array
          uall   = .true.     ! Store all terms in columns
          ulfl   = .true.     ! Assemble upper and lower part
          mtype  = 11         ! PARDISO symmetric indefinite
          if(debug) then
            write(*,2000) '>>>>> PARDISO for non-symmetric ',
     &                    'indefinite matrices activated.'
          endif
          write(iow,2000) '>>>>> PARDISO for non-symmetric ',
     &                    'indefinite matrices activated.'

!         Reset IPARM unsymmetric parameters

          iparm(11) = 1
          iparm(13) = 1

!       Symmetric indefinite solver (DEFAULT)

        else

          ubycol = .false.    ! Sparse storage by columns
          udiag  = .true.     ! Store diagonal array
          uall   = .false.    ! Store all terms in columns
          ulfl   = .false.    ! Assemble upper and lower part
          mtype  = -2         ! PARDISO symmetric indefinite

          if(debug) then
            write(*,2000) '>>>>> PARDISO for symmetric ',
     &                    'indefinite matrices activated.'
          endif
          write(iow,2000) '>>>>> PARDISO for symmetric ',
     &                    'indefinite matrices activated.'

!         Reset IPARM symmetric parameters

          iparm(11) = 0
          iparm(13) = 0

        endif

!       Establish current matrix indices and compute npzs (static variable saved)

        npzsv = npzs
        call iters(npzs,-2)     ! for all row indexed solvers

!       Allocate matrix storage

        if(npzs.ne.npzsv) then
          write(tname,'(4hTANG,i1)') npart
          setvar = palloc(npart,tname,npzs,2)

!         Delete solution array if it exists

          if(np(227).ne.0) then
            setvar = palloc(227,'USOL3', 0, 1)
          endif

!         Release PARDISO memory for new problem

          phase = -1                ! release internal memory
          call pardiso(pt, maxfct, mnum, mtype, phase, neq, ddum, idum,
     1                 idum,idum,nrhs,iparm, istats, ddum, ddum, error)

!         Reset PT pointers

          pt(1:64) = 0

        endif

!       Zero tangent array

        call pzero(hr(np(npart)),npzs)

!       Allocate array to store shifed IA array

        if(np(227).eq.0) then
          setvar = palloc(227,'USOL3', neq+1, 1)
        endif

c       Shifts for IA array

        do i = neq,1,-1
          mr(np(227)+i) = mr(np(225)+i-1) + 1
        end do ! i
        mr(np(227)) = 1

!       Error sets

        error       =  0           ! initialize error flag

      endif

      if(flags(2)) then  ! Factor step

!       Reordering and Symbolic Factorization, This step also allocates
!       all memory that is necessary for factorization

        phase = 11              ! only reordering and symbolic factorization

        call pardiso(pt, maxfct, mnum, mtype, phase, neq,
     &               hr(np(npart)), mr(np(227)), mr(np(226)),
     &               idum, nrhs, iparm, istats, ddum, ddum, error)
        if(debug) then
          write(  *,*) 'Reordering completed ... '
          write(iow,*) 'Reordering completed ... '
          if (error .ne. 0) then
            write(  *,*) 'The following ERROR was detected: ', error
            write(iow,*) 'The following ERROR was detected: ', error
            call plstop()
          endif
          write(  *,*) 'Number of nonzeros in factors = ',iparm(18)
          write(iow,*) 'Number of nonzeros in factors = ',iparm(18)
          write(  *,*) 'Number of factorization MFLOPS = ',iparm(19)
          write(iow,*) 'Number of factorization MFLOPS = ',iparm(19)
        endif

        phase = 22              ! only factorization

        call pardiso(pt, maxfct, mnum, mtype, phase, neq,
     &               hr(np(npart)), mr(np(227)), mr(np(226)),
     &               idum, nrhs, iparm, istats, ddum, ddum, error)

        if(debug) then
          write(  *,*) 'Factorization completed ... '
          write(iow,*) 'Factorization completed ... '
          if (error .ne. 0) then
            write(  *,*) 'The following ERROR was detected: ', error
            write(iow,*) 'The following ERROR was detected: ', error
            call plstop()
          endif
          if((mtype.eq.-2).or.(mtype.eq.11)) then
            write(  *,*) 'Number of zero pivots',iparm(14)
            write(iow,*) 'Number of zero pivots',iparm(14)
          endif
          if(mtype.eq.-2) then
            write(  *,*) 'Number of positive eigenvalues',iparm(22)
            write(iow,*) 'Number of positive eigenvalues',iparm(22)
            write(  *,*) 'Number of negative eigenvalues',iparm(23)
            write(iow,*) 'Number of negative eigenvalues',iparm(23)
          endif
        endif

!       Return error flag to FEAP

        if(error.ne.0) then
          flags(5) = .true.
        endif

      endif

!     Back substitution and iterative refinement

      if(flags(4)) then

        phase = 33              ! only solution step

!       Allocate memory for temporary storage of solution vector

        setvar = palloc(229,'USOL5',neq,2)

        call pardiso(pt, maxfct, mnum, mtype, phase, neq,
     &               hr(np(npart)), mr(np(227)), mr(np(226)),
     &               idum, nrhs, iparm, istats, b, hr(np(229)), error)

        if(debug) then
          write(  *,*) 'Solve completed ... '
          write(iow,*) 'Solve completed ... '
          write(  *,*) 'Number of iterative refinements:',iparm(7)
          write(iow,*) 'Number of iterative refinements:',iparm(7)
        endif

!       Compute energy norm

        aengy = 0.d0
        do i=1,neq
          aengy = aengy +b(i)*hr(np(229)+i-1)
        end do ! i

!       Update solution vector

        do i=1,neq
          b(i) = hr(np(229)+i-1)
        end do ! i

!       Free workspace storage

        setvar = palloc(229,'USOL5',0,2)

      endif

!     Deallocate memory.

      if (flags(5)) then

!       Termination and release of memory

        phase = -1                ! release internal memory
        call pardiso(pt, maxfct, mnum, mtype, phase, neq, ddum, idum,
     1              idum, idum, nrhs, iparm, istats, ddum, ddum, error)

        if(np(225).ne.0) setvar = palloc(225,'USOL1',0,1) ! allocated in uiters
        if(np(226).ne.0) setvar = palloc(226,'USOL2',0,1) ! allocated in uiters
        if(np(227).ne.0) setvar = palloc(227,'USOL3',0,1) ! allocated here

      end if

!     Formats

 2000 format(a,a/a,i4,a)

      end
