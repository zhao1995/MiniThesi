c$Id:$
      subroutine gplas3d(d,eps,epsp,ep,ntm,istrt, sig,dd,dr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]

c     Purpose: Generalized plasticity model for FEAP. Includes
c              elastic plastic isotropic and kinematic hardening

c     Inputs:
c        d      - Array of material constants
c        eps    - Current strains
c        epsp   - Eps_pl(t=t(n)) & back stresses at t(n)
c        ep(*)  - Effective plastic strain & salm at t(n); state
c        ntm    - Number of components
c        istrt  - Start state: 0 = elastic; 1 = last solution

c     Outputs:
c        epsp   - Eps_pl(t=t(n+1)) & back stresses at t(n+1)
c        ep     - Effective platic strain & salm at t(n+1)
c        sig    - Stresses at t(n+1)
c        dd     - Tangent matrix at t(n+1)
c        dr     - Rayleigh damping matrix at t-n+1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'counts.h'
      include  'pconstant.h'
      include  'sdata.h'

      logical   flg_pl, state
      integer   i,j,istrt,ntm
      real*8    alph(6),eps(*),ed(6),epsp(*),sig(*),sigtr(6),saltr(6)
      real*8    d(*),dd(6,6),dr(6,6),en(6),ep(3),sal(6)
      real*8    elam,twog,g1,g2,g3,bulk,g,yield,hh,hiso,hkin,beta
      real*8    p,a1,a2,a3,xk1,a,gam1,gam2,delta,lambda,saltrm,yld
      real*8    theta,salm,salplm,gfac,gfact, dot,yldgpl,dyld

      save

c     Check state for iterations

      if(niter.eq.0) then         ! First iteration in step
        if(istrt.eq.0) then       ! Elastic state requested
          state = .false.
          dyld  =  0.0d0
        else                      ! Last state requested
          state = .false.
          dyld  =  1.0d-08*d(41)
        endif
      else                        ! Not first iteration in step
        state = .true.
        dyld  =  0.0d0
      endif

c     Extract parameters

      g     = d(27)
      bulk  = d(21) - 2.0d0*two3*g
      twog  = g + g

      yield = d(41)
      beta  = d(42) - yield
      delta = d(43)

      hiso  = d(44)
      hkin  = d(45)

      elam  = bulk - twog*one3

      g1    = g + hkin*one3
      g2    = g + (hiso + hkin)*one3
      g3    = g + hiso*one3
      hh    = two3*(hkin + hiso)
      beta  = sqt23*beta
      delta = two3*delta
      gfac  = 1.0d0

c     Extract back stress

      do i = 1,ntm
        alph(i) = epsp(i+ntm)
      end do ! i

c     Compute volumetric and deviatoric strains

      theta = (eps(1) + eps(2) + eps(3))*one3
      do i = 1,3
        ed(i)   = eps(i) - theta
      end do ! i
      do i = 4,min(ntm,6)
        ed(i) = eps(i)*0.5d0
      end do ! i

      theta = theta*3.0d0

c     Compute trial values

      do i = 1,ntm
        sigtr(i) = twog*(ed(i) - epsp(i))
        saltr(i) = sigtr(i) - alph(i)
      end do ! i

c     Compute yield state and trial (s-alpha) norm

      yld    = sqt23 * (yield + hiso*ep(1))
      salm   = ep(2)
      saltrm = sqrt(dot(saltr,saltr,ntm)+dot(saltr(4),saltr(4),ntm-3))

c     1. Trial stress outside yield surface (=> lambda > 0)

      if (saltrm+dyld.gt.yld .and. state) then

        yld = yldgpl(ep(1),salm,saltrm,g,yield,beta,delta,hiso,hkin,
     &               sqt23,lambda)

        if (lambda.gt.0.d0) then
          flg_pl = .true.

c         Update plastic values

          ep(3)  = 1.0d0
          gfact  = 2.0d0 * g1 * lambda
          do i = 1,ntm
            en(i)  = saltr(i)/saltrm
            sal(i) = saltr(i) - gfact*en(i)
          end do ! i

          salplm = saltrm - gfact
          ep(1)  = ep(1)  + sqt23*lambda
          do i = 1,ntm
            alph(i) = alph(i) + two3*hkin*lambda*en(i)
            epsp(i) = epsp(i) + lambda*en(i)
          end do ! i

c         Elasto-plastic tangent matrix: constant preparation

          p    = twog * lambda / saltrm
          a1   = salplm - salm + (delta + two3*hkin)*lambda
          a2   = salplm - yld
          a3   = (delta + hh)*beta
          xk1  = 2.0d0*g2*a1 + (2.0d0*g-delta)*a2 + a3

          a    = twog * (a1 + a2)/xk1
          gam1 = 1.0d0 - p
          gam2 = p - a

c         Elasto-plastic tangent matrix:  matrix generation

          gfac = gam1
          do j=1,ntm
            do i=1,ntm
              dd(i,j) = twog*gam2*en(i)*en(j)
            end do ! i
          end do ! j
        else
          flg_pl= .false.
          ep(3)  = 0.0d0
        end if
      else
        flg_pl= .false.
        ep(3)  = 0.0d0
      end if

c     Update stresses and other parameters depending on type of step

      if (flg_pl) then

        do i=1,ntm
          sig(i)      = sal(i) + alph(i)
          epsp(i+ntm) = alph(i)
        end do ! i
        ep(2)  = salplm

      else

        do i=1,ntm
          sig(i)      = sigtr(i)
          epsp(i+ntm) = alph(i)
        end do ! i
        ep(2)  = saltrm

      end if

c     Add pressure and form final tangents

      g1 = g*gfac
      g2 = g1 + g1
      g3 = g2*one3
      do i = 1,min(3,ntm)
        sig(i) = sig(i) + bulk*theta
        do j = 1,min(3,ntm)
          dd(j,i) = dd(j,i) + bulk - g3
          dr(j,i) = bulk - g*two3
        end do ! j
        dd(i,i) = dd(i,i) + g2
      end do ! i
      do i = 4,min(ntm,6)
        dd(i,i) = dd(i,i) + g1
        dr(i,i) = g
      end do ! i

      end
