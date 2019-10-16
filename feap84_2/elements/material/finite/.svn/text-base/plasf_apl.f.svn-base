cId:$
      subroutine plasf_apl(d,dp,epn,istrt, sig,ee,dd)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    12/12/2012
c       1. Kinematic hardening added                        10/09/1999
c       2. Rewritten / weakly coupled behavior added *LM*   09/07/2012
c                       (i.e. gam, beta_1, j1)
c       3. Add istrt option, remove use of 'tol'            20/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Return mapping algorithm, orthotropic plasticity
c              [Papadopoulos, Lu 2001]

c      Inputs:
c        d(*)    - Parameter array
c        dp(9)   - Papadopoulos-Lu elastic parameters
c        epn(6)  - Plastic strains
c        istrt   - First iteration control

c      Outputs:
c        sig(6)  - Cauchy stress
c        ee(6)   - Elastic strains
c        dd(6,6) - Tangent moduli
c-----[--.----+----.----+----.-----------------------------------------]
c   - YIELD FUNCTION:
c       f = sqrt[ 1/2 k2 (s_22-s_33)^2 + 3/2 k3 (s'_11)^2
c         + 2 k4 (s_12)^2 + 2 k5 (s_23)^2 + 2 k6 (s_13)^2
c         + gam (s_22-s_33) s'_11 ] - sigma

c       with s'_11 = 2/3 s_11 - 1/3 s_22 - 1/3 s_33

c   - POTENTIAL FUNCTION:
c       P = 1/2 K (trE)^2 + 1/2 mu_2 (E_22-E_33)^2 + 3/2 mu_3 (E'_11)^2
c         + 2 mu_4 (E_12)^2 + 2 mu_5 (E_23)^2 + 2 mu_6 (E_13)^2
c         + beta_1 (E_22-E_33) E'_11
c         + beta_5 trE (E_22-E_33) + beta_6 trE E'_11

c       with E'_11 = 2/3 E_11 - 1/3 E_22 - 1/3 E_33
c-----[--.----+----.----+----.-----------------------------------------]
c                                Hill 48
c     Potential Function   |  Yield and Flow Rule   |  Kinematic Hard.
c      dp(1) = K           |   d(51) = F            |   d(151) = h2
c      dp(2) = mu_2        |   d(52) = G            |   d(152) = h3
c      dp(3) = mu_3        |   d(53) = H            |   d(153) = h4
c      dp(7) = mu_4        |   d(54) = L            |   d(154) = h5
c      dp(8) = mu_5        |   d(55) = M            |   d(155) = h6
c      dp(9) = mu_6        |   d(56) = N            |   d(156) = j1
c      dp(4) = beta_1      |                        |
c      dp(5) = beta_5      |                        |
c      dp(6) = beta_6      |                        |
c-----[--.----+----.----+----.-----------------------------------------]
c     Saturation Hardening, d(46)=5    | Swift Power Hardening, d(46)=6
c      d(41) = y0    init yield stress |  d(41) = y0    init yield
c      d(42) = yi    inf  yield stress |  d(42) = K     strength
c      d(43) = beta  saturation exp.   |  d(43) = eps0  initial strain
c      d(44) = H     isotr hardening   |  d(44) = n     exponent
c                                      |  d(45) = H     isotr hardening
c-----[--.----+----.----+----.-----------------------------------------]
c     Strain & Stress
c       epn(i)  i=1,6   plastic strains
c       epn(7)          accumulated plastic strain
c       epn(i)  i=8,13  back stress
c       sig(i)  i=1,6   trial stress tensor compontents
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'counts.h'
      include   'debugs.h'
      include   'eldata.h'
      include   'iofile.h'
      include   'pconstant.h'
      include   'setups.h'

      logical    state, kine
      integer    istrt, i,iter
      real*8     d(*),dp(9),epn(*),sig(6),ee(6),dd(6,6)
      real*8     a3(3),a2(3),n2(3),n3(3),nd(6)
      real*8     lam2,ahat,y0
c     real*8     tol
      real*8     yield,hard,res0,res,hiso,delty,beta,ep, dyld
      real*8     ksw,eps0,nexp
      real*8     s2,s3,s4,s5,s6, s1tr,s2tr,s3tr,s4tr,s5tr,s6tr
      real*8     k2,k3,k4,k5,k6, mu,mu2,mu3,mu4,mu5,mu6
      real*8     mu22,mu32,mu33, mu42, mu52, mu62, beta12
      real*8     j12, k33,k42,k52,k62
      real*8     h2,h3,h4,h5,h6, beta1,beta5,beta6,gam,j1
      real*8     f23(2,2),f4,f5,f6, x23(2,2),x4,x5,x6
      real*8     y23(2,2),y4,y5,y6, temp
      real*8     gs23(2),gs4,gs5,gs6, xs23(2),xs4,xs5,xs6
      real*8     depn(3),ss(3,3),spr(3),rel

      save

c     data       tol  / 1.d-8 /
      data       a2   / 0.0d0, 1.0d0, -1.0d0 /

c     Check state for iterations

      if(niter.eq.0) then         ! First iteration in step
        if(istrt.eq.0) then       ! Elastic state requested
          state = .false.
          dyld  =  0.0d0
        else                      ! Last state requested
          state = .true.
          dyld  =  1.0d-08*d(41)/dp(2)
        endif
      else                        ! Not first iteration in step
        state = .true.
        dyld  =  0.0d0
        if(rank.gt.0) dyld = -1.0d-08*d(41)/dp(2)
      endif

c     Set a3

      a3(1) =  two3
      a3(2) = -one3
      a3(3) = -one3

c     Material Constants

      mu    = dp(2)
      mu2   = 1.0d0                 ! potential function
      mu3   = dp(3)/mu
      mu4   = dp(7)/mu
      mu5   = dp(8)/mu
      mu6   = dp(9)/mu
      beta1 = dp(4)/mu
      beta5 = dp(5)/mu
      beta6 = dp(6)/mu

      beta12= 2.d0*beta1

      mu22  = 2.d0*mu2
      mu32  = 2.d0*mu3
      mu33  = 3.d0*mu3
      mu42  = 2.d0*mu4
      mu52  = 2.d0*mu5
      mu62  = 2.d0*mu6

c     Convert Hill parameters to Papadopoulos-Lu form

      k2    = 0.5d0*(d(52) + d(53)) + 2.d0*d(51)
      k3    = 1.5d0*(d(52) + d(53))
      k4    = d(56)
      k5    = d(54)
      k6    = d(55)
      gam   = 1.5d0*(d(53) - d(53))

      k33   = 3.d0*k3
      k42   = 2.d0*k4
      k52   = 2.d0*k5
      k62   = 2.d0*k6

      h2    = d(151)/mu             ! kinematic hardening moduli
      h3    = d(152)/mu
      h4    = d(153)/mu
      h5    = d(154)/mu
      h6    = d(155)/mu
      j1    = d(156)/mu

      j12   = 2.d0*j1

      kine  = nint(d(150)).eq.1     ! Indicator for kinematic hardening

c     Hardening Law. Determine current value of yield stress

      y0    = d(41)/mu
      if(nint(d(46)).eq.5) then     ! Saturation hardening
        delty = (d(42) - d(41))/mu
        if(delty.eq.0.0d0) then
          beta = 0.0d0
        else
          beta = d(43)
        endif
        hiso  =  d(44)/mu
        yield = y0 + hiso*epn(7) + delty*(1.0d0-exp(-beta*epn(7)))
      elseif(nint(d(46)).eq.6) then ! Swift power law hardening
        ksw   = d(42)/mu
        eps0  = d(43)
        nexp  = d(44)
        hiso  = d(45)/mu
        yield = hiso*epn(7) + ksw*(eps0 + epn(7))**nexp
      endif

c     Scale stress tensor

      sig = sig/mu

c     Trial Stress Factors (from stgm6.f)

      s1tr =(sig(1) + sig(2) + sig(3))*one3
      s2tr = sig(2) - sig(3)
      s3tr = sig(1) - s1tr
      s4tr = sig(4)               ! tensorial (not engineering) shear
      s5tr = sig(5)
      s6tr = sig(6)

      if(kine) then
        s1tr = s1tr - (epn(8) + epn(9) + epn(10))*one3
        s2tr = s2tr - (epn(9) - epn(10))
        s3tr = s3tr - (two3*epn(8) - (epn(9) + epn(10))*one3)
        s4tr = s4tr - epn(11)
        s5tr = s5tr - epn(12)
        s6tr = s6tr - epn(13)
      end if

c     Check for yield

      res0  = ( 0.5d0*(k2*s2tr*s2tr + k33*s3tr*s3tr)
     &      + ( k42*s4tr*s4tr + k52*s5tr*s5tr + k62*s6tr*s6tr)
     &      +   gam*s2tr*s3tr)/(yield**2) - 1.0d0


c     Plastic return mapping, local iteration for consistency parameter

c     if(res0.gt.0.0d0) then
      if(res0.gt.dyld .and. state) then

        if(kine) then
          y23(1,1) = (mu22 + h2)     *k2  + two3*(beta1 +j1)  *gam
          y23(1,2) = (mu22 + h2)     *gam +      (beta12+j1)  *k3
          y23(2,1) = two3*(beta1+j1) *k2  + (two3*mu3+one3*h3)*gam
          y23(2,2) = two3*(beta1+j1) *gam +      (mu32+h3)    *k3
          y4       = (mu42 + h4)     *k4
          y5       = (mu52 + h5)     *k5
          y6       = (mu62 + h6)     *k6
        else
          y23(1,1) = two3*beta1 *gam + mu22      *k2
          y23(1,2) = mu22       *gam + beta12    *k3
          y23(2,1) = two3*mu3   *gam + two3*beta1*k2
          y23(2,2) = two3*beta1 *gam + mu32      *k3
          y4       = mu42 *k4
          y5       = mu52 *k5
          y6       = mu62 *k6
        end if

c       Iteration loop start

        res  = res0
        lam2 = 0.0d0
        iter = 0
 10     iter = iter + 1

          ep    = epn(7) + sqt23*lam2
          yield = y0 + hiso*ep + delty*(1.0d0-exp(-beta*ep))

          if(nint(d(46)).eq.5) then     ! Saturation hardening yield

            yield = y0   + hiso*ep + delty*(1.0d0-exp(-beta*ep))
            hard  = sqt23*(hiso    + delty*  beta*exp(-beta*ep))

          elseif(nint(d(46)).eq.6) then ! Swift power law hardening

            yield =        hiso*ep +      ksw*(eps0 + ep)**nexp
            hard  = sqt23*(hiso    + nexp*ksw*(eps0 + ep)**(nexp-1.0d0))

          endif

          f23(1,1) = yield + lam2*y23(1,1)
          f23(1,2) =         lam2*y23(1,2)
          f23(2,1) =         lam2*y23(2,1)
          f23(2,2) = yield + lam2*y23(2,2)
          call invert(f23,2,2)

          f4 = 1.0d0/(yield + lam2*y4)
          f5 = 1.0d0/(yield + lam2*y5)
          f6 = 1.0d0/(yield + lam2*y6)

          s2 = f23(1,1)*s2tr + f23(1,2)*s3tr
          s3 = f23(2,1)*s2tr + f23(2,2)*s3tr
          s4 = f4*s4tr
          s5 = f5*s5tr
          s6 = f6*s6tr

          res = 0.5d0*(k2*s2*s2 + k33*s3*s3)
     &        + (k42*s4*s4 + k52*s5*s5 + k62*s6*s6)
     &        +  gam*s2*s3 - 1.0d0

          x23(1,1) = 0.5d0*( gam*f23(2,1) + k2 *f23(1,1) )
          x23(1,2) = 0.5d0*( gam*f23(2,2) + k2 *f23(1,2) )
          x23(2,1) = 0.5d0*( gam*f23(1,1) + k33*f23(2,1) ) ! symmetric
          x23(2,2) = 0.5d0*( gam*f23(1,2) + k33*f23(2,2) )
          x4       = k42*f4
          x5       = k52*f5
          x6       = k62*f6

          xs23(1)  = x23(1,1)*s2 + x23(1,2)*s3
          xs23(2)  = x23(2,1)*s2 + x23(2,2)*s3
          xs4      = x4*s4
          xs5      = x5*s5
          xs6      = x6*s6

          ahat = s2*((hard+y23(1,1))*xs23(1) +       y23(2,1) *xs23(2))
     &         + s3*(      y23(1,2) *xs23(1) + (hard+y23(2,2))*xs23(2))
     &         + s4*(hard+y4)*xs4 + s5*(hard+y5)*xs5 + s6*(hard+y6)*xs6

          lam2 = lam2 + res/(ahat*2.0d0)

          if(debug.and.n.eq.1) then
            rel = res/res0
            write(iow,'(1x,a,i4,2e14.4)')
     &          'Iteration, Residual, Relative:',iter,res,rel
            if(ior.lt.0)
     &        write(*,'(1x,a,i4,2e14.4)')
     &          'Iteration, Residual, Relative:',iter,res,rel
          end if

c         if(abs(res).gt.res0*1.0d-12.and.abs(res).gt.1.0d-14) go to 10
          if(abs(res).gt.res0*1.0d-10.and.abs(res).gt.1.0d-10) go to 10

c       Interation loop end

c       Update plastic and total strain

        gs23(1) = 0.5d0*lam2*(gam*s3 + k2 *s2)
        gs23(2) = 0.5d0*lam2*(gam*s2 + k33*s3)
        gs4     = lam2*k4*s4
        gs5     = lam2*k5*s5
        gs6     = lam2*k6*s6

        depn(1:3) = gs23(1)*a2(:) + gs23(2)*a3(:)
        epn(1:3)  = epn(1:3) + depn(1:3)
        ee(1:3)   = ee(1:3)  - depn(1:3)
        epn(4)    = epn(4) + gs4
        epn(5)    = epn(5) + gs5
        epn(6)    = epn(6) + gs6
        ee(4)     = ee(4)  - gs4
        ee(5)     = ee(5)  - gs5
        ee(6)     = ee(6)  - gs6
c       epn(7)    = epn(7) + sqt23*lam2*(1.0d0-tol)
        epn(7)    = epn(7) + sqt23*lam2

c       Update kinematic hardening (back-stress tensor) (note: scaled with mu)

        if (kine) then
          epn(8:10) = epn(8:10) +  (h2*gs23(1) + two3*j1*gs23(2)) *a2(:)
     &                          + (j12*gs23(1) +      h3*gs23(2)) *a3(:)
          epn(11)   = epn(11) + h4*gs4
          epn(12)   = epn(12) + h5*gs5
          epn(13)   = epn(13) + h6*gs6
        end if ! back-stress tensor

c       Update Stress Tensor

        n2(:)    =   mu2*a2(:) + beta1*a3(:) + beta5
        n3(:)    = beta1*a2(:) +  mu33*a3(:) + beta6
        sig(1:3) = sig(1:3) - 2.0d0*gs23(1)*n2(:) - two3*gs23(2)*n3(:)
        sig(4)   = sig(4) - mu42*gs4
        sig(5)   = sig(5) - mu52*gs5
        sig(6)   = sig(6) - mu62*gs6

c       Consistent Tangent Operator

        do i = 1, 3
          dd(i,1:3) = dd(i,1:3) - mu*lam2*four3* (
     &          3.0d0*x23(1,1)*n2(i)*n2(:) +      x23(1,2)*n2(i)*n3(:)
     &              + x23(2,1)*n3(i)*n2(:) + one3*x23(2,2)*n3(i)*n3(:) )
        end do ! i
        dd(4,4) = dd(4,4) - mu*lam2*x4*mu4*mu4
        dd(5,5) = dd(5,5) - mu*lam2*x5*mu5*mu5
        dd(6,6) = dd(6,6) - mu*lam2*x6*mu6*mu6

c       Norm-One Modification

        nd(1:3) = 2.0d0*xs23(1)*n2(:) + two3*xs23(2)*n3(:)
        nd(4)   = xs4*mu4
        nd(5)   = xs5*mu5
        nd(6)   = xs6*mu6

        temp = mu*(yield-hard*lam2)/ahat
        do i = 1, 6
          dd(:,i) = dd(:,i) - nd(:)*nd(i)*temp
        end do ! i

      end if ! plastic return mapping

c     Scale stresses

      sig = sig*mu

c     Compute Cockcraft-Lantham Parameter

      ss(1,1) = sig(1)
      ss(2,2) = sig(2)
      ss(3,3) = sig(3)
      ss(1,2) = sig(4)
      ss(2,1) = sig(4)
      ss(2,3) = sig(5)
      ss(3,2) = sig(5)
      ss(3,1) = sig(6)
      ss(1,3) = sig(6)
      call eig3(ss,spr,i)

      epn(14) = epn(14) + max(spr(1),spr(2),spr(3),0.d0)*sqt23*lam2

      end
