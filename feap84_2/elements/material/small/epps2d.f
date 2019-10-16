c$Id:$
      subroutine epps2d(d,eps,epsp,alp,epn,istrt, sig,dd,dr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c     2. Correct start state for inelastic begin            19/03/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Elastic plastic plane stress kinematic/isotropic
c                hardening model: Return map solution

c      Inputs:
c         d(*)      - Array of material constants
c         eps(*)    - Total strain
c         epsp(*)   - Plastic strain at t-n
c         alp(*)    - Back-stress at t-n
c         epn(*)    - Equivalent plastic strain at t-n and state
c         istrt     - Start state: 0 = elastic; 1 = last solution

c      Outputs:
c         sig(*)    - Stresses at t-n+1
c         epsp(*)   - Plastic strain at t-n+1
c         alp(*)    - Back-stresses at t-n+1
c         epn(*)    - Equivalent plastic strain at t-n+1
c         dd(*,*)   - Material tangent matrix at t-n+1
c         dr(*,*)   - Rayleigh damping matrix at t-n+1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'counts.h'
      include  'pconstant.h'

      logical   state
      integer   i,j,istrt
      real*8    yldtr,rad2,xlam0,fact,a1,a2, xkprim,hprim,beta, dyld
      real*8    eps(*),sig(*),sigtr(3),alp(*),d(*),dd(6,6),dr(6,6)
      real*8    en(4),eta(3),epsp(*),epn(2),dum(6), hard2d,ylds2d

      save

c     Check state for iterations

      if(niter.eq.0) then         ! First iteration in step
        if(istrt.eq.0) then       ! Elastic state requested
          state = .false.
          dyld  = 0.0d0
        else                      ! Last state requested
          state = .true.
          dyld  = 1.e-08*d(41)
        endif
      else                        ! Not first iteration in step
        state = .true.
        dyld  = 0.0d0
      endif

c     Set-up elastic moduli

      call dmat2d(d,0.d0,dd,dum)
      do i = 1,4
        do j = 1,4
          dr(j,i) = dd(j,i)
        end do ! j
      end do ! i

c     Compute trial stress

      sigtr(1) = dd(1,1)*(eps(1)-epsp(1)) + dd(1,2)*(eps(2)-epsp(2))
      sigtr(2) = dd(2,1)*(eps(1)-epsp(1)) + dd(2,2)*(eps(2)-epsp(2))
      sigtr(3) = dd(4,4)*(eps(4)-epsp(3))

c     Compute yield state

      eta(1)  = sigtr(1) - alp(1)
      eta(2)  = sigtr(2) - alp(2)
      eta(3)  = sigtr(3) - alp(3)
      yldtr   = (eta(1)**2+eta(2)**2-eta(1)*eta(2))*one3 + eta(3)**2
      rad2    = hard2d(d,epn(1),xkprim,hprim)
      xlam0   = 0.0d0

c     Compute scale factor and compute elastic-plastic tangent

      if (yldtr+dyld.gt.rad2 .and.state) then

        epn(2) = 1.0d0
        dm     = ylds2d(d,sigtr,alp,eta,dd,beta,epn(1),xlam0,hprim)

c       Compute plastic strains

        a1      = (2.d0*eta(1) - eta(2))*one3
        a2      = (2.d0*eta(2) - eta(1))*one3
        epsp(1) = epsp(1) + xlam0*a1
        epsp(2) = epsp(2) + xlam0*a2
        epsp(3) = epsp(3) + xlam0*2.d0*eta(3)

c       Compute a "normal" and finish computation of tangent

        en(1) = dd(1,1)*a1 + dd(1,2)*a2
        en(2) = dd(2,1)*a1 + dd(2,2)*a2
        en(3) = 0.0d0
        en(4) = dd(4,4)*eta(3)*2.0d0
        fact  = 1.d0/max(abs(en(1)),abs(en(2)),abs(en(4)))
        en(1) = en(1)*fact
        en(2) = en(2)*fact
        en(4) = en(4)*fact
        beta  = beta *fact
        fact  = sqrt(en(1)*a1 + en(2)*a2 + en(4)*eta(3)*2.d0 + beta)
        if(fact.le.0.0d0) write(*,*) ' ERROR IN FACT',fact
        do i = 1,4
          en(i) = en(i)/fact
        end do ! i
        do i = 1,4
          do j = 1,4
            dd(i,j) = dd(i,j) - en(i)*en(j)
          end do ! j
        end do ! i

      else
        epn(2) = 0.0d0
      endif

c     Set stresses

      sig(1) = sigtr(1)
      sig(2) = sigtr(2)
      sig(3) = 0.0d0
      sig(4) = sigtr(3)

      end
