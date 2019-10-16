c$Id:$
      subroutine plas1d(d,ta, eps,h1,istrt, sig,dd)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Take alpha from d(47) instead of d(3)            10/07/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 1-d plasticity constitutive equation

c      Inputs:
c         d(*)   - Parameters
c         ta     - Temperature
c         eps    - Strain
c         hl(*)  - History parameters
c         istrt  - Start state: 0 = elastic; 1 = last solution

c      Outputs:
c         sig    - Stress
c         dd(*)  - Modul
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'counts.h'

      logical   state, start
      integer   istrt, nseg
      real*8    ta,eps, sig,dd(2), d(*),h1(*)

      save

c     Check state for iterations

      if(niter.eq.0) then         ! First iteration in step
        if(istrt.eq.0) then       ! Elastic state requested
          state = .false.
          start = .false.
        else                      ! Last state requested
          state = .false.
          if(h1(3).gt.0.0d0) then
            start = .true.
          else
            start = .false.
          endif
        endif
      else                        ! Not first iteration in step
        state = .true.
        start = .false.
      endif

      nseg = nint(d(130))
      if(nseg.eq.0) then
        call pllh1d(d,ta, eps,h1,start, sig,dd, state)
      else
        call plsh1d(d,d(131),nseg, ta, eps,h1,start, sig,dd, state)
      endif

      end

      subroutine pllh1d(d,ta, eps,h1,start, sig,dd, state)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  1-d linear kinematic, saturation isotropic hardening

c      Inputs:
c         d(41)    - Initial yield stress
c         d(42)    - Yield at infinite accumulated plastic strain
c         d(43)    - Exponential hardening constant
c         d(44)    - Linear isotropic hardening modulus
c         d(45)    - Linear Kinematic hardening modulus
c         ta       - Temperature
c         eps      - Strain at t_n+1
c         h1(1)    - Plastic strain at t_n               (epn)
c         h1(2)    - Accumulated plastic strain at t_n   (epp)
c         sig      - Initial stress specified
c         start    - Start request

c      Outputs:
c         sig      - Stress at t_n+1
c         dd(1)    - Elastic-Plastic modulus at t_n+1
c         h1(1)    - Plastic strain at t_n+1             (epn)
c         h1(2)    - Accumulated plastic strain at t_n+1 (epp)
c         h1(3)    - State indicator: Elastic = 0; Plastic = 1

c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'iofile.h'

      logical  noconv, state,start
      integer  count
      real*8   ta,eps, sig,dd(2), epp,epn,alp, ff,yld,etol,d(*),h1(*)
      real*8   epstr,sigtr,cc, dyld,Rs,rf,nn,det,dsig,dlam,lambda

      data     etol /1.d-08/

c     Extract effective plastic strain

      epn   = h1(1)
      epp   = h1(2)

c     Compute trial stress

      sigtr = sig + d(21)*(eps - epn - d(47)*ta)

c     Check yield

      alp   = d(45)*epn
      ff    = abs(sigtr-alp)
      yld   = d(42) + (d(41) - d(42))*exp(-d(43)*epp) + d(44)*epp

c     Plastic state

      if((ff.gt.yld .and. state) .or. start) then

        lambda =  0.0d0
        sig    =  sigtr

        cc     =  1.d0/d(21)
        epstr  =  cc*sigtr

c       Newton solution for state

        noconv = .true.
        count  =  0

        do while(noconv)

c         Compute Newton parameters

          count =  count + 1
          dyld  = (d(42) - d(41))*exp(-d(43)*epp)
          yld   =  d(42) - dyld + d(44)*epp
          dyld  =  d(43)*dyld + d(44) + d(45)
          ff    = abs(sig - alp)
          nn    = (sig - alp)/ff
          Rs    =  epstr - cc*sig - lambda*nn
          rf    =  yld - ff
          det   =  1.d0/(dyld*cc + 1.d0)

c         Compute increments

          dsig  =  det*(dyld*Rs + nn*rf)
          dlam  =  det*(  nn*Rs - cc*rf)

c         Update variables

          lambda = lambda + dlam
          sig    = sig    + dsig

          epn    = h1(1)  + lambda*nn
          epp    = h1(2)  + lambda
          alp    = d(45)*epn

c         Check convergence

          if(abs(d(41)-d(42)).le.etol*abs(d(41))) then ! Linear hard.
            noconv = .false.
          elseif(abs(dlam).lt.etol*abs(lambda)) then   ! Newton converge
            if(abs(dsig).lt.etol*abs(sig))  then
              noconv = .false.
            endif
          elseif(count.gt.25) then                     ! Count limit
            if(.not.start) then
              if(ior.lt.0) then
                write(*,2000) lambda,dlam
              endif
              write(iow,2000) lambda,dlam
            endif
            noconv = .false.
          endif
        end do ! while noconv

c       Elasto-plastic modulus

        dd(1) = d(21)*(1.d0 - det)

c       Update history variables

        h1(1) = epn
        h1(2) = epp
        h1(3) = 1.d0

c     Elastic state

      else

c       Elastic modulus

        sig   = sigtr
        dd(1) = d(21)
        h1(3) = 0.d0

      endif

      dd(2)  = d(21)

2000  format(' *WARNING* No convergence in PLAS1D'/
     &       '           lambda =',1p,1e12.4,' dlambda =',1p,1e12.4)
      end

      subroutine plsh1d(d,table,nseg, ta, eps,h1,start, sig,dd, state)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Plasticity: 1-d table lookup hardening/yield

c     Input parameters
c        d(*)       - Material parameters
c        table(3,*) - Table Values: 1 = plastic strain
c                                   2 = yield stress
c                                   3 = hardening parameter
c        nseg       - Number of segments in table
c        ta         - Temperature
c        eps        - Strain
c        h1(*)      - History variables
c        start      - Start flag
c     Output parameters
c        sig        - Stress
c        dd         - Elasto-plastic modulus
c        state      - State indicator
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'

      logical    noconv, state,start
      integer    nseg,count,i
      real*8     ta,eps, sig,dd(2), epp, epn, alp,d(*),table(3,*),h1(*)
      real*8     ee,cc, rsig,ralp,rf,nn,det,dsig,dalp,dlam,lambda
      real*8     a11,a12,a13,a22,a23, ff,yld, iso, kin, etol

      save

      data       etol /1.d-08/

c     Extract effective plastic strain

      epn = h1(1)
      epp = h1(2)
      alp = h1(4)

c     Compute trial stress

      sig = sig + d(21)*(eps - epn - d(47)*ta)

      i = 1
      call tyld1d(nseg,table, epp, yld, iso, kin, i)

      ff = abs(sig - alp)
      nn = (sig - alp)/ff

c     Plastic state

      if((ff.gt.yld .and. state) .or. start) then

        lambda =  0.0d0
        ee     =  d(21)
        cc     =  1.d0/ee

c       Newton solution for state

        noconv = .true.
        count  =  0

        do while(noconv)

c         Compute Newton parameters

          count =  count + 1

c         Compute parameters

          call tyld1d(nseg,table, epp, yld, iso, kin, i)

          ff   =  abs(sig - alp)
c         nn   = (sig - alp)/ff
          rsig =  eps - h1(1) - cc*sig - lambda*nn
          if(kin.gt.0.0d0) then
            ralp = -(alp - h1(4))/kin + lambda*nn
          else
            ralp = 0.0d0
          endif
          rf   = -ff + yld
          det  =  1.d0/(ee + iso + kin)
          a11  =  ee*(iso + kin)
          a12  =  ee*kin
          a13  =  ee*nn
          a22  =  kin*(ee + iso)
          a23  = -kin*nn

c         Compute increments

          dsig = (a11*rsig + a12*ralp + a13*rf)*det
          dalp = (a12*rsig + a22*ralp + a23*rf)*det
          dlam = (a13*rsig + a23*ralp -1.d0*rf)*det

c         Update variables

          lambda = lambda + dlam
          sig    = sig    + dsig
          alp    = alp    + dalp

          epn    = h1(1)  + lambda*nn
          epp    = h1(2)  + lambda

c         Check convergence

          if(  abs(dlam).le.etol*abs(lambda)) then
            if(abs(dsig).le.etol*abs(sig) .and.
     &         abs(dalp).le.etol*abs(alp))  then
              noconv = .false.
            endif
          endif

c         Check count limit

          if(count.gt.25) then
            if(.not.start) then
              if(ior.lt.0) then
                write(*,*) ' *WARNING* No convergence in PLSH1D',dlam,
     &                     dsig,dalp
              endif
              write(iow,*) ' *WARNING* No convergence in PLSH1D',dlam,
     &                     dsig,dalp
            endif
            noconv = .false.
          endif

        end do ! while noconv

c       Elasto-plastic modulus

        dd(1) = a11*det

c       Update history variables

        h1(1) = epn
        h1(2) = epp
        h1(3) = 1.d0

c     Elastic state

      else

c       Elastic modulus

        dd(1) = d(21)
        h1(3) = 0.d0

      endif
      dd(2) = d(21)

      end

      subroutine tyld1d(nseg,table, epp, yld, iso, kin, i)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Table Look-up for Yield and hardening moduli

c     Input parameters
c        nseg       - Number of segments in table
c        table(3,*) - Table Values: 1 = plastic strain
c                                   2 = yield stress
c                                   3 = hardening parameter
c        epp        - Accumulated plastic strain
c     Output parameters
c        yld        - Yield stress
c        iso        - Isotropic hardening
c        kin        - Kinematic hardening
c        i          - Table location
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'

      logical    flag
      integer    nseg,i
      real*8     epp, del, yld, iso,kin, table(3,*)

c     Find table entries which bound the current epp

      flag = .true.
      i    = 1
      do while (flag)
        if(epp.lt.table(1,i)) then  ! epp < entry (i-1 < epp < i)
          flag = .false.
        elseif(i.ge.nseg) then      ! epp > last entry
          flag = .false.
          i    = nseg
        else
          i    = i + 1
        endif
      end do ! while

c     Error: epp < first entry

      if(i.eq.1) then
        write(iow,4000) epp,table(1,1)
        if(ior.lt.0) then
          write(*,4000) epp,table(1,1)
        endif
        call plstop()
      endif

c     Interpolate table

      del =  1.d0/(table(1,i) - table(1,i-1))
      iso = (table(2,i) - table(2,i-1))*del
      yld = (table(2,i-1)*table(1,i) - table(2,i)*table(1,i-1))*del
     &    +  iso*epp
      kin = (table(3,i-1)*table(1,i) - table(3,i)*table(1,i-1))*del
     &    + (table(3,i) - table(3,i-1))*del *epp

c     Format

4000  format(/' *ERROR: PLAS1D: epp =',1p,1e12.5,' < table(1) = ',
     &        1p,1e12.5)

      end
