c$Id:$
      subroutine dparam(ct,lct)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c        1. Add start option for Newmark                    02/01/2007
c        2. Revise setting of 'nrk', 'nrc' for 'gen1'       06/07/2007
c           algorithm
c        3. Add standard HHT option: trans hht alpha        23/07/2008
c           and Euler implicit     : trans eule
c        4. Restrict write of theta to 3 in format 4000     30/11/2009
c        5. Allow start with 'a_n+1' = 'a_n'                14/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set integration parameters for Time-stepping Algorithms

c      Inputs:
c        ct(3)   := Algorithmic parameters
c        lct     := 'off ' (noi=0)  Standard Static
c                   'newm' (noi=1)  Classical Newmark
c                   'back' (noi=2)  Backward  Euler
c                   'alph' (noi=3)  Conserving HHT (0.7 < alpha < 1)
c                   'hht ' (noi=3)  Standard   HHT (0.7 < alpha < 1)
c                   'expl' (noi=4)  Explicit Newmark
c                   'cons' (noi=5)  Conserving Newmark (Momentum&Energy)
c                   'stat' (noi=6)  Generalized Mid-point Static
c                   'gen1' (noi=7)  Generalized Mid-point First Order
c                   'cent' (noi=8)  Central Difference Explicit
c                   'bdf2' (noi=9)  Backward Difference Formula (2)
c                   'eule' (noi=10) Euler implicit

c                   'user' (noi=-1) User time integration routine
c                   'init'          Initialize the nrt variable

c      Outputs:
c        ct(3)   := (possibly) redefined algorithmic parameters
c        theta(i):= Set identical to ct(i), i=1,2,3.
c        noi     := Number of the specific integrator (see above)
c        nrm     := Acceleration vector (pointer) in urate(nneq,*)
c        nrk     := Solution vector (pointer)
c        nrc     := Velocity vector (pointer)
c        nrt     := Maximun number of vectors in urate(nnneq,*)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'counts.h'
      include  'ddata.h'
      include  'debugs.h'
      include  'dyndat.h'
      include  'fdata.h'
      include  'gltran.h'
      include  'iofile.h'
      include  'part0.h'
      include  'part7.h'
      include  'sdata.h'
      include  'tdata.h'

      logical   pcomp
      character lct*4, u0(4)*1, u1(4)*3
      integer   i, ord, ntot(10)
      real*8    ct(3)

      save

      data      u0 / 'u'  , 'v'  , 'a'  , 'a'   /
      data      u1 / 'u_n', 'v_n', '0.0', 'a_n' /

c     Set maximum number of vectors for each integration type
      data ntot/   2   ,   1 ,    4 ,    2    ,   5   ,   3   ,  2,
c               Newmark, Back,   HHT, Explicit, Conser, Static, First
c                       Euler  Alpha            ving            Order

     &             2   ,   3   ,   2  /
c               Central Backward Euler
c            Difference Differ.  Impl.

c     Initialize the nrt variable

      if(pcomp(lct,'init',4)) then

        nrt = max(ntot(1),ntot(2),ntot(3),ntot(4),ntot(5),ntot(6),
     &            ntot(7),ntot(8),ntot(9),ntot(10))
        i   = 0
        call uparam(ct,nrk,nrc,nrm,ord,i,0)
        nrt = max(nrt,i) + 2

      else

        nrk = 0
        nrc = 0
        nrm = 1

        dynflg = .true.

c       Classical Newmark-beta method

        if(pcomp(lct,'    ',4).or.pcomp(lct,'newm',4)) then

c         ct(1) = beta  ;  ct(2) = gamma

          noi = 1
          nrk = 0
          nrc = 1
          nrm = 2
          ord = 2

          if(ct(1).eq.0.0d0) ct(1) = 0.25d0
          if(ct(2).eq.0.0d0) ct(2) = 0.5d0
          idyn0 = max(1,min(4,nint(ct(3))))

          write(iow,2001) ct(1),ct(2),u0(idyn0),u1(idyn0)
          if(ior.lt.0) then
            write(*,2001) ct(1),ct(2),u0(idyn0),u1(idyn0)
          endif
          ct(3) = 1.0d0

c       Backward Euler for first order ode (e.g. heat transfer)

        elseif(pcomp(lct,'back',4)) then

c        'nrk' is solution, 'nrc' is first rate term
c        'nrm' is set to same address as 'nrc'
c         No theta(i) values are required (i.e., theta(1) is 1.0!)

          noi = 2
          nrk = 0
          nrc = 1
          nrm = 1
          ord = 1

          ct(1) = 1.0d0
          ct(2) = 0.0d0
          ct(3) = 1.0d0
          write(iow,2002)
          if(ior.lt.0) write(*,2002)

c       Conserving HHT alpha-method

        elseif(pcomp(lct,'alph',4)) then

c         ct(1) = beta  ;  ct(2) = gamma ; ct(3) = alpha

          noi = 3
          nrk = 3
          nrc = 4
          nrm = 2
          ord = 2

          if(ct(1).eq.0.0d0) ct(1) = 0.5d0
          if(ct(2).eq.0.0d0) ct(2) = 1.0d0
          if(ct(3).eq.0.0d0) ct(3) = 0.5d0

          write(iow,2003) ct(1),ct(2),ct(3)
          if(ior.lt.0) write(*,2003) ct(1),ct(2),ct(3)

c       Standard HHT alpha-method

        elseif(pcomp(lct,'hht',3)) then

c         ct(1) = alpha

          noi = 3
          nrk = 3
          nrc = 4
          nrm = 2
          ord = 2

          if(ct(1).eq.0.0d0) ct(1) = 1.0d0
          ct(3) = ct(1)
          ct(1) = 0.25d0*(2.d0 - ct(3))**2
          ct(2) = 1.50d0 - ct(3)

          write(iow,2003) ct(1),ct(2),ct(3)
          if(ior.lt.0) write(*,2003) ct(1),ct(2),ct(3)

c       EXPLICIT Newmark method (beta=0)

        elseif(pcomp(lct,'expl',4)) then

c         ct(1) = gamma  or;  ct(2) = gamma

          noi = 4
          nrk = 0
          nrc = 1
          nrm = 2
          ord = 2

          if(ct(1).eq.0.0d0) ct(1) = 0.5d0
          if(ct(2).eq.0.0d0) ct(2) = ct(1)
          ct(1) = 0.0d0
          ct(3) = 1.0d0

          write(iow,2004) ct(2)
          if(ior.lt.0) write(*,2004) ct(2)

c       MOMENTUM & ENERGY conserving algorithm

        elseif(pcomp(lct,'cons',4)) then

c         ct(1) = beta  ;  ct(2) = gamma ; ct(3) = alpha

          noi = 5
          nrk = 3
          nrc = 4
          nrm = 5
          ord = 2

          if(ct(1).eq.0.0d0) ct(1) = 0.5d0
          if(ct(2).eq.0.0d0) ct(2) = 1.0d0
          if(ct(3).eq.0.0d0) ct(3) = 0.5d0

          write(iow,2005) ct(1),ct(2),ct(3)
          if(ior.lt.0) write(*,2005) ct(1),ct(2),ct(3)

c       Static Solution Procedure for Generalized Mid-point Algorithms

        elseif(pcomp(lct,'stat',4)) then

c         ct(1) = alpha  ;  ct(2) = None ; ct(3) = alpha

          noi = 6
          nrk = 3
          nrc = 4
c         nrm = 5
          nrm = 2
          ord = 0

          if(ct(1).eq.0.0d0) ct(1) = 0.5d0
          if(ct(3).eq.0.0d0) ct(3) = ct(1)
          ct(1) = 1.0d0
          ct(2) = 1.0d0

          write(iow,2006) ct(3)
          if(ior.lt.0) write(*,2006) ct(3)

          dynflg = .false.

c       First Order Solution: Generalized Mid-point Algorithm

        elseif(pcomp(lct,'gen1',4)) then

c         ct(1) = alpha  ;  ct(2) = None ; ct(3) = alpha

          noi = 7
          nrk = 3
          nrc = 1
          nrm = 2
          ord = 1

          if(ct(1).eq.0.0d0) ct(1) = 0.5d0
          if(ct(3).eq.0.0d0) ct(3) = ct(1)
          ct(1) = 1.0d0
          ct(2) = 1.0d0

          write(iow,2007) ct(3)
          if(ior.lt.0) write(*,2007) ct(3)

c       Central Difference scheme for second order equations

        elseif(pcomp(lct,'cent',4)) then

          noi = 8
          nrk = 0
          nrc = 1
          nrm = 2
          ord = 2

          ct(1) = 1.0d0
          ct(2) = 1.0d0
          ct(3) = 1.0d0

          write(iow,2008)
          if(ior.lt.0) write(*,2008)

c       First Order Solution: Backward Difference Formula (2)

        elseif(pcomp(lct,'bdf2',4)) then

c         ct(1) = alpha  ;  ct(2) = None ; ct(3) = alpha

          noi = 9
          nrk = 0
          nrc = 1
          nrm = 1
          ord = 1

          ct(1) = 1.0d0
          ct(2) = 1.0d0

          write(iow,2009)
          if(ior.lt.0) write(*,2009)

          nstepa = nstep

c       Backward Euler for 2nd order problems

        elseif(pcomp(lct,'eule',4)) then

          noi = 10
          nrk = 0
          nrc = 1
          nrm = 2
          ord = 2

          write(iow,2010)
          if(ior.lt.0) write(*,2010)

c       User time integration routine

        elseif(pcomp(lct,'user',4)) then

          noi    = -1
          call uparam(ct,nrk,nrc,nrm,ord,i,1)
          dynflg        = i.gt.0
          fl(9)         = dynflg
          flp(9,npart)  = dynflg

c         Set order for default ones

          do i = 1,ndf
            if(npart.eq.ndfp(i)) then
              ndfo(i) = min(ord,ndog(i))
            endif
          end do ! i

c       Standard static algorithm (noi = 0)

        elseif(pcomp(lct,'off',3)) then

          noi     =  0
          nrk     =  0
          nrc     =  0
          nrm     =  0
          ord     =  0

          ct(1)   = 0.0d0
          ct(2)   = 0.0d0
          ct(3)   = 1.0d0

          dynflg       = .false.
          fl(9)        = .false.
          flp(9,npart) = .false.

          gtan(1) = 1.0d0
          gtan(2) = 0.0d0
          gtan(3) = 0.0d0

          write(iow,2000)
          if(ior.lt.0) write(*,2000)

c       ERROR - write message

        else

          if(ior.lt.0) then
            write(*,3000) lct
          else
            write(iow,3000) lct
            write(ilg,3000) lct
            call plstop()
          endif
        endif

c       Transfer values to 'theta' in common /ddata/

        do i = 1,4
          theta(i) = ct(i)
        end do ! i
        alpha  =  ct(3)
        numint =  noi

c       Set order for default ones

        do i = 1,ndf
          if(npart.eq.ndfp(i)) then
            ndfo(i) = min(ord,ndog(i))
          endif
        end do ! i

c       Debug information

        if(debug) then
          write(iow,4000) noi,nrk,nrc,nrm,ord,(theta(i),i=1,3)
          if(ior.lt.0) then
            write(*,4000) noi,nrk,nrc,nrm,ord,(theta(i),i=1,3)
          endif
        endif

      endif

2000  format(/' Standard Static Algorithm'/)

2001  format(/' Newmark Parameters:',
     &        ' Beta = ',f9.4,' Gamma = ',f9.4/
     &    20x,' Start with ',a1,'_n+1 = ',a3/1x)

2002  format(/' Backward Euler for First Order Systems.'/1x)

2003  format(/' Conserving HHT Parameters'/
     &        '  Beta = ',f9.4,' Gamma = ',f9.4,' Alpha = ',f9.4/1x)

2004  format(/' Newmark Explicit:  Gamma = ',f9.4/1x)

2005  format(/' Conserving Algorithm Parameters'/
     &        '  Beta = ',f9.4,' Gamma = ',f9.4,' Alpha = ',f9.4/1x)

2006  format(/' Static Generalized Mid-Point Configuration:',
     &        ' Alpha = ',f9.4/1x)

2007  format(/' First Order Generalized Mid-Point Integration:',
     &        ' Alpha = ',f9.4/1x)

2008  format(/' Central Difference Explicit Integration:')

2009  format(/' Backward Difference Formula (BDF2) for First Order',
     &        ' Systems.'/1x)

2010  format(/' Backward Euler for 2nd Order Equations')

3000  format(/' *ERROR* DPARAM: ',a,' Not an implemented method'/1x)

4000  format(/' Dynamic Solution Variable Locations'/
     &       5x,'Solution method    (noi) = ',i3/
     &       5x,'Stiffness location (nrk) = ',i3/
     &       5x,'Damping   location (nrc) = ',i3/
     &       5x,'Mass      location (nrm) = ',i3/
     &       5x,'ODE order          (ord) = ',i3/
     &       5x,'Theta(1) factor          = ',1p,1e12.5/
     &       5x,'Theta(2) factor          = ',1p,1e12.5/
     &       5x,'Theta(3) factor          = ',1p,1e12.5)
      end
