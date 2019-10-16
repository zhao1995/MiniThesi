c$Id:$
      subroutine mises(d,eps,epsp,epp,ntm,istrt, sig,dd,dr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Mises (J2) plasticity with isotropic and kinematic hardening

      implicit  none

      include  'counts.h'
      include  'elauto.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'sdata.h'
      include  'setups.h'
      include  'tdata.h'

      logical   conv,state
      integer   i, j,istrt, ntm, count, mm,m1
      real*8    d(*),eps(*),epsp(*),epp(*), sig(*), dd(6,6),dr(6,6)
      real*8    alp(6), ep(6), en(6), xi(6)
      real*8    aa,bb,cc, press, theta, dlam, lam, xin, tolc, dot
      real*8    k, g, gbar, twog, hi,hi23,his23,hk,hk23, epss, dyld
      real*8    r0,rinf,rn,rb, beta, expb,expl, bbar, s23
      real*8    ff,phi,dphi, Tbar, vbar

      save

      data      tolc / 1.d-10 /

c     Check state for iterations

      if(niter.eq.0) then         ! First iteration in step
        if(istrt.eq.0) then       ! Elastic state requested
          state = .false.
          dyld  =  0.0d0
        else                      ! Last state requested
          state = .true.
          dyld  =  1.0d-08*d(41)
        endif
      else                        ! Not first iteration in step
        state = .true.
        dyld  =  0.0d0
        if(rank.gt.0) dyld = -1.0d-08*d(41)
      endif

c     Set parameters

      s23   = sqrt(two3)

      g     = d(27)
      k     = d(21) - 2.d0*two3*g

      hi    = d(44)
      hk    = d(45)

      r0    = s23*d(41)
      rinf  = s23*d(42)
      beta  =     d(43)

c     Viscoplastic parameters

      if(dt.gt.0.0d0) then
        Tbar  = d(180)/dt
      else
        Tbar  = 0.0d0
      endif
      mm = nint(d(181))
      m1 = mm - 1

c     Compute constants

      hk23  = two3*hk
      hi23  = two3*hi
      his23 =  s23*hi

      twog  = 2.d0*g

      expb  = exp(-beta*epp(1))
      bbar  =  s23*beta*(rinf - r0)*expb
      gbar  = twog + hi23 + hk23

c     Compute volumetric and deviatoric strains

      theta = (eps(1) + eps(2) + eps(3))*one3
      do i = 1,3
        ep(i)   = eps(i) - theta
      end do ! i
      do i = 4,ntm
        ep(i) = eps(i)*0.5d0
      end do ! i

      theta = theta*3.0d0

c     Compute trial values

      do i = 1,ntm
        sig(i) = sig(i) + twog*(ep(i) - epsp(i))
        alp(i) = hk23*epsp(i)
        xi(i)  = sig(i) - alp(i)
      end do ! i

c     Compute trial norm of stress - back stress

      xin = sqrt(dot(xi,xi,ntm) + dot(xi(4),xi(4),ntm-3))
      rn  = rinf + his23*epp(1)
      rb  = (r0  - rinf)*expb

c     Check yield

      if(xin+dyld .gt. (rn+rb) .and. state) then

c       Compute strain factor

        if(r0.gt.0.0d0) then
          epss   = r0/twog
        else
          epss   = 1.e-4
        endif

c       Compute viscoplastic or plasticity consistency

        conv   = .false.
        count  =  0

c       Plasticity

        if(d(180).eq.0.0d0)  then
          vbar = 0.0d0
          lam  = (xin - rn - rb) / (gbar + bbar)
          do while (.not.conv .and. count.le.25)
            expl = exp(-s23*beta*lam)
            dlam = (xin - rn - rb*expl - gbar*lam)/(gbar + bbar*expl)
            lam  = lam + dlam
            if(abs(r0-rinf).le.tolc*r0) then
              conv = .true.
            elseif(abs(dlam) .lt. tolc*abs(lam) .or.
     &         abs(lam)  .lt. tolc*epss)          then
              conv = .true.
            endif
            count = count + 1
          end do ! while

c       Viscoplasticity

        else
          lam = 0.0d0
c         lam = (xin - rn - rb) / (gbar + bbar)
          do while (.not.conv .and. count.le.25)
            expl = exp(-s23*beta*lam)
            ff   = xin - rn - rb*expl - gbar*lam
            if(m1.gt.0) then
              dphi  = (ff/r0)**m1 / r0
            else
              dphi  = 1.d0/r0
            endif
            phi  = dphi*ff
            dphi = dphi*mm

            dlam = (phi - Tbar*lam)/(dphi*(gbar + bbar*expl) + Tbar)
            lam  = lam + dlam
            if(abs(dlam) .lt. tolc*abs(lam) .or.
     &         abs(lam)  .lt. tolc*epss)          then
              conv = .true.
            endif
            count = count + 1
          end do ! while

c         Mask lambda to be positive

          lam  = max(0.0d0,lam)
          vbar = Tbar/dphi

        endif

c       Warning: Not converged

        if(.not.conv .and. niter.gt.0) then
          write(  *,*) '  *WARNING* No convergence in MISES'
          write(iow,*) '  *WARNING* No convergence in MISES'
          write(iow,*) '   lam ',lam,' dlam ',dlam,' count ',count
c         call plstop()
        endif

c       Set auto time stepping factor

        rmeas = max(rmeas,2.5d0*rvalu(1)*s23*lam/epss)

c       Compute normal vector

        do i = 1,ntm
          en(i) = xi(i)/xin
        end do ! i

        if(lam.gt.0.0d0) then
          bb     = twog*lam/xin
          aa     = twog*(1.d0 - bb)
          bb     = twog*(bb - twog/(gbar + bbar*expl + vbar))
          epp(2) = 1.0d0
        else
          aa     = twog
          bb     = 0.0d0
          epp(2) = 0.0d0
        endif

c       Compute plastic tangent from normal

        do i = 1,ntm
          cc = bb*en(i)
          do j = 1,ntm
            dd(i,j) = cc*en(j)
          end do ! j
        end do ! i

c     Set for elastic increment

      else

        epp(2) = 0.0d0
        do i = 1,ntm
          en(i) = 0.0d0
        end do ! i

        lam = 0.0d0
        aa  = twog

        do i = 1,ntm
          do j = 1,ntm
            dd(i,j) = 0.0d0
          end do ! j
        end do ! i

      endif

c     Compute deviator stress, plastic strain, & accumul plastic strain

      do i = 1,ntm
        sig(i)  = sig(i)  - twog*lam*en(i)
        epsp(i) = epsp(i) + lam*en(i)
      end do ! i
      epp(1) = epp(1) + s23*lam

c     Add pressure

      press = k*theta
      do i = 1,3
        sig(i) = sig(i) + press
      end do ! i

c     Compute tangent moduli

      cc = k - aa*one3
      bb = k - twog*one3

      do i = 1,3
        do j = 1,3
          dd(j,i) = dd(j,i) + cc
          dr(j,i) =           bb
        end do ! j
        dd(i  ,i  ) = dd(i  ,i  ) + aa
        dr(i  ,i  ) = dr(i  ,i  ) + twog
      end do ! i
      do i = 4,ntm
        dd(i,i) = dd(i,i) + 0.5d0*aa
        dr(i,i) = g
      end do ! i

      end
