c$Id:$
      subroutine umacr10(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Authors : D.S. Bindel, S. Govindjee                    21/06/2005
c     Modifies: R.L. Taylor
c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct the names of deleted USER files          31/08/2007
c       2. Delete redundant deletion arrays                 03/09/2007
c       3. For 'symm' option allow computation of all the   07/09/2007
c          eigenvalues/vectors.
c       4. Move setting of fp(1) outside test for lump      12/09/2007
c       5. Add checks on number of mass terms and set the   13/09/2007
c          mode 1 options.  Report zero/negative mass term.
c       6. Change to umacr10                                25/01/2012
c       7. Add set of 'vneq' and 'pfeapb.h'                 16/02/2013
c          Change 'mmode' to 'emode'
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Serial ARPACK interface to FEAP

c      Inputs:
c         lct        - Command character parameters
c                      'symm'etric: (K - lambda M)*V = 0
c         ctl(3)     - Command numerical parameters
c                      1 = Number eigenpairs to compute
c                      2 = Maximum number of iterations
c                      3 = Tolerance for eigenvalues

c      Outputs:
c         hr(np(76)) - Eigenvalues
c         hr(np(77)) - Eigenvectors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'debugs.h'
      include  'evdata.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'rdata.h'
      include  'umac1.h'
      include  'p_int.h'
      include  'p_point.h'
      include  'part0.h'
      include  'pfeapb.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   pcomp, setval,palloc
      character lct*15
      real*8    ctl(3),stol

c     ARPACK PARAMETERS BEGIN

      integer   maxitr, emode, nn

c     ARPACK PARAMETERS END

      save

c     Set command word to: ARPAck

      if(pcomp(uct,'ma10',4)) then      ! Usual    form
        uct = 'arpa'                    ! Specify 'arpa'ck
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else

c       Set number of eigenvalues and size

        mf     = min(neq,max(1,nint(ctl(1))))
        mq     = min(neq,max(20,2*mf))
        maxitr = max(300,nint(ctl(2)))
        vneq   = neq

c       Set convergence tolerance

        if(ctl(3).ne.0.0d0) then
          stol = ctl(3)
        else
          stol = max(1.d-12,tol)
        endif

c       Symmetric eigensolution

        if(pcomp(lct,'symm',4) .or. pcomp(lct,'lump',4)
     &                         .or. pcomp(lct,'    ',4)) then

c         Count number of diagonal mass entries

          if(fl(1)) then
            fp(1) = np(npart+8)
          else
            fp(1) = np(npart+12)
          endif
          call numass(hr(fp(1)),neq,mq)

c         Reduce number of eigenvalues and output numbers

          mf = min(mf,mq)

          if(ior.lt.0) then
            write(*,2000) mf,maxitr,stol
          endif

c         Lumped mass algorithm not requiring matrix factors

          if(pcomp(lct,'lump',4)) then

c           Check for correct number of mass entries

            if(mq.lt.neq) then
              write(iow,3000) mq,neq
              write(ilg,3000) mq,neq
              if(ior.lt.0) then
                write(*,3000) mq,neq
                return
              endif
              call plstop()
            endif

c           Set algorithm to mode = 1

            emode = 1
            setval = palloc(157,'USER7',neq,2)    ! Sqrt M
            point  = np(157)
            call masssqr(hr(point),hr(np(12+npart)), neq)

c         Lumped/Consistent mass form requires factors of stiffness

          else
            point = np(1)
            emode = 3
          endif

c         Reduce requested number: Arpack will compute anyway

          if(mf.eq.mq) then
            mf  = mf - 1
            nn  = 1
          else
            nn  = 0
          endif

c         Diagnostic prints in debug mode

          if(ior.lt.0 .and. debug) then
            write(  *,2001) 'Symmetric',mf+nn,mq,shift,maxitr,neq
          elseif(debug) then
            write(iow,2001) 'Symmetric',mf+nn,mq,shift,maxitr,neq
          endif

c         Allocate eigenpair and working arrays for symmetric problem

          setval = palloc( 76,'EVAL ', 2*mq        , 2)  ! d
          setval = palloc( 77,'EVEC ', mq*neq      , 2)  ! v
          setval = palloc(158,'USER8', neq         , 2)  ! resid
          setval = palloc(159,'USER9', mq*(mq+8)   , 2)  ! workl
          setval = palloc(160,'USER0', 3*neq       , 2)  ! workd

          call arfeaps(emode,      hr(np( 76)), hr(np( 77)),
     &                 hr(point)  ,hr(np(158)), hr(np(159)),
     &                 hr(np(160)),mf, mq, maxitr, shift,stol)
          if(emode.eq.1) then
            setval = palloc(157,'USER7',  0,2)    ! Delete mass sqrt
          endif

c       Error mode

        else
           if(ior.lt.0) then
             write(  *,3001)
             return
           else
             write(iow,3001)
             write(ilg,3001)
             call plstop()
          endif
        endif

c       Delete scratch arrays

        setval = palloc(158,'USER8', 0,2)  ! resid
        setval = palloc(159,'USER9', 0,2)  ! workl
        setval = palloc(160,'USER0', 0,2)  ! workd

      endif

c     Format

2000  format(/'  Computing',i4,' Eigenvalues: Max iter =',i5,
     &        ' Tol =',1p,1e12.4/1x)


2001  format(/10x,a,' Eigenproblem',/,
     &        12x,'No. vecs requested  ',i6,/
     &        12x,'No. vecs being used ',i6,/
     &        12x,'Shift:              ',1p,1e15.8,/
     &        12x,'Maximum # iterations',i6,/
     &        12x,'Number eqns. (neq)  ',i6)

3000  format(5x,'** ERROR ** Mass matrix has',i8,' terms for',i8,
     &          'equations')

3001  format(5x,'** ERROR in ARPACK ** Specify Option: SYMM or ITER')

      end

      subroutine masssqr(msqrt,mass, neq)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Compute reciprocal of square root of mass

c     Inputs:
c       mass(*)  - Lumped mass
c       neq      - Number of equations

c     Outputs:
c       msqrt(*) - Reciprocal of square root of mass
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'iofile.h'

      logical    ermass
      integer    neq,n
      real*8     msqrt(*), mass(*)

      ermass = .false.
      do n = 1,neq
        if(mass(n).gt.0.0d0) then
          msqrt(n) = 1.d0/sqrt(mass(n))
        else
          if(ermass) then
            write(iow,3001) n,mass(n)
            write(ilg,3001) n,mass(n)
            if(ior.lt.0) then
              write(*,3001) n,mass(n)
            endif
          else
            write(iow,3000) n,mass(n)
            write(ilg,3000) n,mass(n)
            if(ior.lt.0) then
              write(*,3000) n,mass(n)
            endif
          endif
          ermass = .true.
        endif
      end do ! n

c     Exit if error

      if(ermass) then
          call plstop()
      endif

c     Format

3000  format(5x,'** ERROR ** Lumped mass has zero or negative entries.'/
     &       12x,'Term:',i8,' Value:',1p,1e15.5)

3001  format(12x,'Term:',i8,' Value:',1p,1e15.5)

      end
