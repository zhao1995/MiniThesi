c$Id:$
      subroutine umacr4(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Authors:  D.S. Bindel, S. Govindjee                    21/06/2005
c     Modified: R.L. Taylor                                  14/12/2005
c       Original version                                     01/11/2006
c       1. Delete 'prt' from argument, add 'print.h'         21/04/2011
c       2. Changed from 'prt' to 'prnt' for pring flg        05/01/2013
c       3. Remove pneq just use numpeq                       07/01/2013
c       4. Add 'parpacka.h'; call to count number of masses  11/04/2013
c          Set vector to unity; add format 2000
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Parallel ARPACK interface to FEAP using PETSc.
c                instructions.
c      Notes:
c         (1) Solution with consistent mass requires solution of
c             equations: K * z = x.
c         (2) Use of multi-grid solver requires equations to be blocked
c             for more efficient solution.
c         (3) Lumped mass can use either blocked or unblocked form of
c             stiffness (K).  If frequecies are above ~1000Hz may need
c             to modify routine pminvsqr.F

c      Inputs:
c         lct = lump - Command character parameters
c                      Mode 1: Symmetric with lumped mass
c                        (M^-1/2*K*M^-1/2 - lambda I) * X = 0
c                         V = M^-1/2 * X
c         lcr = xxxx   Mode 3: Symmetric with consistent mass
c                        (K - lambda M) * V = 0
c                      N.B. xxxx can be anything except 'lump'.
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
      include  'evdata.h'
      include  'iofile.h'
      include  'pfeapb.h'
      include  'pointer.h'
      include  'print.h'
      include  'rdata.h'
      include  'sdata.h'
      include  'umac1.h'

      include  'comblk.h'

      include  'parpacka.h'

      logical   pcomp, setval,palloc
      character lct*15
      real*8    ctl(3),stol

c     ARPACK STUFF BEGIN

      integer   maxitr, n

c     ARPACK STUFF END

      save

c     Set command word

      if(pcomp(uct,'mac4',4)) then      ! Usual    form
        uct = 'parp'                    ! Specify 'parp'ack
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else

        mf     = min(numteq,max(1,nint(ctl(1))))
        mq     = min(numpeq,max(20,2*mf))
        maxitr = max(300,nint(ctl(2)))
        if(ctl(3).ne.0.0d0) then
          stol = ctl(3)
        else
          stol = max(1.d-12,tol)
        endif

        maxitr = 1

        vneq   = numnp*ndf

c       Check for maximum number of eigenvalues

        setval = palloc(157,'USER7', numpeq*2,2)
        do n = 0,numpeq-1
          hr(np(157)+n) = 1.0d0
        end do ! n

        parmassfl = .true.
        call paropm(numpeq,hr(np(157)),hr(np(157)+numpeq))
        parmassfl = .false.

        setval = palloc(157,'USER7', 0,2)
        mf     = min(mf,nummasst)  ! Reduce number for fewer masses
        write(iow,2000) mf
        if(ior.lt.0) then
          write(*,2000) mf
        endif

c       Symmetric eigensolve

        if(pcomp(lct,'lump',4)) then             ! Lumped mass form
          mmode = 1
          setval = palloc(157,'USER7', numpeq, 2)  ! Temporary
          call pminvsqr(hr(np(157)),mr(np(245)))
          setval = palloc(157,'USER7',    0, 2)  ! Temporary
        else                                     ! Consistent mass form
          mmode = 3
        endif

        if(ior.lt.0) then
          write(*,2001) 'Symmetric',mf,mq,shift,maxitr,numpeq,mmode
        endif
        write(iow,2001) 'Symmetric',mf,mq,shift,maxitr,numpeq,mmode

        setval = palloc( 76,'EVAL ', 2*mq,       2)  ! d - eigenvalues
        setval = palloc( 77,'EVEC ', mq*vneq  ,  2)  ! v - eigenvectors
        setval = palloc(158,'USER8', vneq,       2)  ! resid
        setval = palloc(159,'USER9', mq*(mq+8),  2)  ! workl
        setval = palloc(160,'USER0', 3*vneq,     2)  ! workd

        call parfeaps(mmode,      hr(np( 76)), hr(np( 77)), hr(np(158)),
     &               hr(np(159)), hr(np(160)), numpeq,  vneq, mf,  mq,
     &               maxitr,      shift,       stol,  prnt)

c       Delete scratch arrays

        setval = palloc(160,'USER0', 0,2)  ! workd
        setval = palloc(159,'USER9', 0,2)  ! workl
        setval = palloc(158,'USER8', 0,2)  ! resid

      end if

c     Format

2000  format(/10x,'No. vecs reduced to  ',i6)

2001  format(/10x,a,' Eigenproblem',/,
     &        12x,'No. vecs requested   ',i6,/
     &        12x,'No. vecs being used  ',i6,/
     &        12x,'Shift:      ',1p,1e15.8,/
     &        12x,'Maximum # iterations ',i6,/
     &        12x,'Number eqns. (numpeq)  ',i6:/
     &        12x,'ARPACK solution mode ',i6)

      end
