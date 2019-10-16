c$Id:$
      subroutine plasfd(d, detf, f, epp, be, epl,
     &                  ntm, istrt, dd, sig, isw, state)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c     2. Add YY to argument of yfunct, store epp and yield  06/02/2007
c     3. Remove local define of 'One2' (in pconstant.h)
c     4. Introduce 'tolb' for checking equal roots          15/02/2008
c        Set tolb = 1.d-08
c     5. Evolve effective plastic strain by sqrt(2/3) gam   10/08/2009
c     6. Add 'gam' to yfunct argument list                  06/05/2010
c     7. Add 'epp(2)' to store the SIGMA_n norm of gplast   18/05/2010
c        Add 'ssn' to call of yfunct.
c     8. Add df to argument list                            06/12/2011
c        Use f(3,3,4) to compute solution values.
c     9. Convert to use displacement gradient instead of F  07/12/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Finite Deformation Isotropic I1-J2-J3 Plasticity Models
c              in Principal Logarithmic Stretches
c              Incremental  Lagrangian  Formulation

c     Inputs:
c       d(*)      -  Material parameters
c       detf      -  Jacobian determinant (F)        at t_n+1
c       f(3,3,*)  -  Deformation gradient            at t_n+1
c       istrt     -  Start state: 0 = elastic; 1 = last solution

c                 History Variables
c       epp(1)    -  Cumulative plastic strain       at t_n
c       epp(2)    -  | SIGMA_n | general plasticity  at t_n
c       be(*)     -  Left Cauchy-Green tensor        at t_n
c       epl(*)    -  Plastic Strain for Hardening    at t_n

c     Outputs:
c       sig(*)    -  Cauchy stress tensor
c       dd(6,6)   -  Cauchy (spatial) elastic moduli

c                 History Variables
c       epp(1)    -  Cumulative plastic strain       at t_n+1
c       epp(2)    -  | SIGMA_n+1 | general plastic   at t_n+1
c       be(*)     -  Left  Cauchy-Green tensor       at t_n+1
c       epl(*)    -  Plastic Strain for Hardening    at t_n+1
c-----[--.----+----.----+----.-----------------------------------------]
c     MATERIAL CONSTANTS:
c       d(21)    -  Elastic Bulk  modulus
c       d(22)    -  Elastic Shear modulus

c       d(44)    -  Isotropic Hardening Modulus - H_iso [x sqrt(2/3)]
c       d(45)    -  Kinematic Hardening Modulus - H_kin [x sqrt(2/3)]

c       d(46)    -  1 = VON MISES      [ sqrt(2*J2)]
c       d(41)    -  Yield stress (ep = 0)(Y0)           [x sqrt(2/3)]
c       d(42)    -  Yield stress (ep =00)(Yi)           [x sqrt(2/3)]
c       d(43)    -  Delta value  (ep = 0)(delta)

c       d(46)    -  2 = DRUCKER PRAGER [ sqrt(2*J2) + (1/3) * aa * I1]
c       d(41)    -  Yield stress         (Y0)
c       d(42)    -  aa parameter         (aa)

c       d(46)    -  3 = PRAGER LODE [sqrt(2*J2) + sqrt(27/2)*bb*J3/J2 ]
c       d(41)    -  Yield stress         (Y0)
c       d(42)    -  Lode angle parameter (bb)

c       d(46)    -  4 = Generalized plasticity
c       d(41)    -  Yield stress (ep = 0)(Y0)           [x sqrt(2/3)]
c       d(42)    -  Yield stress (ep =00)(Yi)           [x sqrt(2/3)]
c       d(43)    -  Delta value  (ep = 0)(delta)

c     VARIABLES:
c       vol       - volumetric strain   -   vol = ( eps : 1 ) / 3
c       be(*)     - Left Cauchy-Green tensor
c       nn_t(3,3) - principal directions (by columns) tensor form
c       ll2(3)    - squares of principal stretches =  lambda^2

c       tau(3)    - principal values   total    Kirchhoff stress
c       tt(3)     - principal values deviatoric Kirchhoff stress
c       pp        - pressure         volumetric Kirchhoff stress

c       dtde(3,3) - Kirchhoff stress derivative
c                   dtde(a,b)   = [d tau_a/d (lambda_b^2) ]*lambda_b^2
c                               = [d tau_a/d ( eps_b ) ]

c       dd_l(6,6)   - spatial elastic moduli in principal basis

c     FUNCTIONS:
c       yield     - yield function value
c       yfunct    - yield function
c                    flg = .false.  only yield function
c                    flg = .true.   yield function and derivatives
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'counts.h'
      include 'iofile.h'
      include 'pconstant.h'
      include 'rdata.h'
      include 'setups.h'

      logical  conv, state
      integer  ntm,istrt, a, b, rot, it, itmax, isw, p1(6), p2(6)
      real*8   d(*), be(*),epl(*),dd(6,6), sig(*)
      real*8   f(3,3,*),fi(3,3),detfr,detf, detr
      real*8   pp, kk , gg, TwoG, K3inv, G2inv, Ypr, YY, dyld
      real*8   ll2_tr(3), ll2(3), tau(3), tt(3)
      real*8   eps_tr(3), th_tr, vol_tr, ee_tr(3), alp_n(3)
      real*8   eps_e(3),         vol_e,  ee_e(3) , alp(3), ta(3)
      real*8   Y0, I1, J2, J3, f1, f2, f3, f11, f22, f33, f12, f13, f23
      real*8   gam, Hk, Hkr, aatr, epse, xe
      real*8   err, dsol(7),res(7), tres(7,7),fss(3,3)
      real*8   dtde(3,3), ff(6,6), dd_l(6,6), epp(2), eppn, ssn
      real*8   nn_t(3,3), nn(3), aa(3), bb(3)
      real*8   tolb,tolc, d_el, xx,f112,f123, yield,yfunct

      real*8   delf(3,3), delb(3,3), dfcp(3,3)

      data     p1    /1,2,3,1,2,3/
      data     p2    /1,2,3,2,3,1/
      data     dd_l  / 36*0.0d0 /
      data     tolb  / 1.0d-08 /

c     INITIALIZE HISTORY VARIABLES: ISW = 14

      if (isw.eq.14) then        ! Initialize non-zero quantities
        return
      end if

c     READING MATERIAL PARAMETERS  &  CONSTANT SETUP

      kk    = d(21)
      gg    = d(22)
      Y0    = d(41)  ! N.B. Radius of yield = sqrt(2/3) sig_y

      Hk    = d(45)  ! N.B. Kinematic hardening = 2/3 * H_kin
      if(Hk.gt.0.0d0) then
        Hkr = 1.d0/Hk
      else
        Hkr = 0.0d0
      endif

      TwoG  = 2.0d0 * gg

      itmax = 50                      ! Iteration maximum #
      tolc  = 1.d-09                  ! Tolerance for convergence

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
        if(rank.gt.0) dyld  = -1.0d-08*d(41)
      endif

c     KINEMATIC COMPUTATIONS  (Incremental Lagrangian formulation)

c     Inverse of f(3,3,2) * J_n

      fi(1,1) = f(2,2,2)*f(3,3,2) - f(2,3,2)*f(3,2,2)
      fi(2,1) = f(2,3,2)*f(3,1,2) - f(2,1,2)*f(3,3,2)
      fi(3,1) = f(2,1,2)*f(3,2,2) - f(2,2,2)*f(3,1,2)

      fi(1,2) = f(3,2,2)*f(1,3,2) - f(3,3,2)*f(1,2,2)
      fi(2,2) = f(3,3,2)*f(1,1,2) - f(3,1,2)*f(1,3,2)
      fi(3,2) = f(3,1,2)*f(1,2,2) - f(3,2,2)*f(1,1,2)

      fi(1,3) = f(1,2,2)*f(2,3,2) - f(1,3,2)*f(2,2,2)
      fi(2,3) = f(1,3,2)*f(2,1,2) - f(1,1,2)*f(2,3,2)
      fi(3,3) = f(1,1,2)*f(2,2,2) - f(1,2,2)*f(2,1,2)

      detfr   = 1.d0/(f(1,1,2)*fi(1,1)
     &              + f(1,2,2)*fi(2,1)
     &              + f(1,3,2)*fi(3,1))

c     Compute (G_n+1 - G_n)/F_n^-1

      do b = 1,3
        do a = 1,3
          delf(a,b) =((f(a,1,3) - f(a,1,4)) * fi(1,b)
     &              + (f(a,2,3) - f(a,2,4)) * fi(2,b)
     &              + (f(a,3,3) - f(a,3,4)) * fi(3,b))*detfr
        end do ! a
      end do ! b

      delb(1,1) = be(1)
      delb(2,2) = be(2)
      delb(3,3) = be(3)

      delb(1,2) = be(4)
      delb(2,1) = be(4)

      delb(2,3) = be(5)
      delb(3,2) = be(5)

      delb(3,1) = be(6)
      delb(1,3) = be(6)

c     Compute Elastic left Cauchy-Green tensor:    be = F * Cp * F^T

      do b = 1,3
        do a = 1,3
          dfcp(a,b) = delf(a,1) * delb(1,b)
     &              + delf(a,2) * delb(2,b)
     &              + delf(a,3) * delb(3,b)
        end do ! a
      end do ! b
      do b = 1,3
        do a = 1,3
          nn_t(a,b) = delf(a,b) + delf(b,a) + delb(a,b)
     &              + dfcp(a,b) + dfcp(b,a)
     &              + delf(a,1) * delf(b,1)
     &              + delf(a,2) * delf(b,2)
     &              + delf(a,3) * delf(b,3)
     &              + dfcp(a,1) * delf(b,1)
     &              + dfcp(a,2) * delf(b,2)
     &              + dfcp(a,3) * delf(b,3)
        end do ! a
      end do ! b

c     Save 'be' for elastic case

      be(1) = nn_t(1,1)
      be(2) = nn_t(2,2)
      be(3) = nn_t(3,3)
      be(4) = nn_t(1,2)
      be(5) = nn_t(2,3)
      be(6) = nn_t(3,1)

c     Compute principal stretches and directions

      call eig3(nn_t,ll2_tr,rot)

c     COMPUTE TRIAL KIRCHHOFF STRESS  ( pressure and deviator )

      do a = 1, 3
        if(abs(ll2_tr(a)).gt.0.001) then
          ee_tr(a) = 0.5d0*log( 1.0d0 + ll2_tr(a) ) ! log(lambda(a)^TR)
        else
          xe       = ll2_tr(a)
          ee_tr(a) = 0.5d0*(xe - one2*xe**2
     &                         + one3*xe**3
     &                         - one4*xe**4
     &                         + one5*xe**5
     &                         - one6*xe**6)
        endif
      end do ! a
      th_tr = ee_tr(1) + ee_tr(2) + ee_tr(3)
      do a = 1,3
        ee_tr(a) = ee_tr(a) - one3*th_tr   ! Trial deviators
      end do ! a

      pp   = kk * th_tr                    ! Pressure: K*th_tr
      eppn = epp(1)
      ssn  = epp(2)                        ! Norm of (ss - alp)_n

      do a = 1, 3
        tt(a)    = TwoG   * ee_tr(a)       ! Trial deviatoric stress
        tau(a)   = tt(a)  + pp
        alp(a)   = Hk * epl(a)
        alp_n(a) = alp(a)
      end do ! a

c     Deviatoric: ta =  tt - alp_dev

      aatr = (alp(1) + alp(2) + alp(3))*one3
      do a = 1,3
        ta(a)   = tt(a) - alp(a) + aatr    ! Trial SIGMA = ss - alp
      end do ! a

c     CHECK ELASTIC / PLASTIC STEP

c     Compute stress invariant and yield function

      I1 =  3.0d0 * (pp - aatr)
      J2 = ( ta(1)**2 + ta(2)**2 + ta(3)**2 ) * one2
      J3 = ( ta(1)**3 + ta(2)**3 + ta(3)**3 ) * one3

      gam   = 0.0d0
      yield = yfunct(d,epp,I1,J2,J3, f1,f2,f3,
     &               f11,f22,f33,f12,f13,f23,Ypr,YY,gam,ssn,.false.)
c    &               f11,f22,f33,f12,f13,f23,Ypr,YY,.false.)

c     Check yield

      if ( (yield+dyld .gt. 0.0d0) .and. state ) then

c     PLASTIC STEP  -->  RETURN MAP

c       Update solution if gamma greater than zero from yfunct

        if(gam.gt.0.0d0 .and. nint(d(46)).eq.4) then

          f2 = 1.0d0/(sqrt(ta(1)**2 +ta(2)**2 +ta(3)**2))
          do a = 1, 3
            nn(a)  = f2 * ta(a)
            tt(a)  = tt(a)  - TwoG * gam * nn(a)     ! Deviatoric stress
            tau(a) = tt(a)  + pp
            alp(a) = alp(a) + Hk   * gam * nn(a)
          end do ! a

c         Deviatoric: ta =  tt - alp_dev

          aatr = (alp(1) + alp(2) + alp(3))*one3
          do a = 1,3
            ta(a)   = tt(a) - alp(a) + aatr    ! Trial SIGMA = ss - alp
          end do ! a

c         Compute stress invariant and yield function

          I1 =  3.0d0 * (pp - aatr)
          J2 = ( ta(1)**2 + ta(2)**2 + ta(3)**2 ) * one2
          J3 = ( ta(1)**3 + ta(2)**3 + ta(3)**3 ) * one3
        endif

        conv   = .false.
        it     =  1

        K3inv  = one3 /  kk
        G2inv  = one2 /  gg

        vol_tr = th_tr * one3
        vol_e  = K3inv * pp

        do a = 1, 3
          eps_tr(a) = ee_tr(a) + vol_tr
          ee_e(a)   = G2inv * tt(a)
          eps_e(a)  = ee_e(a)  + vol_e
        end do ! a
        epse = 1.d0/(abs(eps_tr(1)) + abs(eps_tr(2)) + abs(eps_tr(3)))

        d_el = ( K3inv - G2inv ) * one3

        do while ((.not.conv).and.(it.le.itmax))

          yield = yfunct(d,epp,I1,J2,J3, f1,f2,f3,
     &                   f11,f22,f33,f12,f13,f23,Ypr,YY,gam,ssn,.true.)
c    &                   f11,f22,f33,f12,f13,f23,Ypr,YY,.true.)

          xx  = f1 - f3 * two3 * J2

          do a = 1, 3
            nn(a)    = xx + ( f2 + f3 * ta(a) ) * ta(a)
            res(a)   = eps_e(a) - eps_tr(a) + gam * nn(a)
            res(a+3) = (alp_n(a) - alp(a))*Hkr + gam * nn(a)
          end do ! a

          res(7) = yield

          err  = (abs(res(1)) + abs(res(2)) + abs(res(3)))*epse
     &         +  abs(res(7))/Y0
          conv = err .lt. tolc

c         Construct local tangent matrix

          do a = 1, 3
            aa(a) = f23   * ta(a) + f13
            bb(a) = ta(a) * ta(a) - two3 * J2
          end do ! a

          f112 = f11 - one3 * f2
          f123 = f12 - two3 * f3

          do a = 1, 3
            do b = a, 3
              fss(a,b) = ( f112 + f22 * ta(a) * ta(b)
     &                 +   f33 * bb(a) * bb(b)
     &                 +  f123 * ( ta(b) + ta(a) )
     &                 + aa(a) * bb(b) + bb(a) * aa(b) ) * gam
              fss(b,a) = fss(a,b)
            end do ! b
            fss(a,a) = fss(a,a) + gam * ( f2 + 2.0d0 * f3 * ta(a) )
          end do ! a

          do a = 1, 3
            do b = 1, 3
              tres(a,b) = fss(a,b) + d_el
            end do ! b
            tres(a,a) = tres(a,a) + G2inv
          end do ! a

c         Kinematic and Isotropic Hardening

          if(Hk.gt.0.0d0) then

            do a = 1,3
              do b = 1,3
                tres(a,b+3)   = -fss(a,b)
                tres(a+3,b)   = -fss(a,b)
                tres(a+3,b+3) =  fss(a,b)
              end do ! b
c             tres(a+3,a+3) = tres(a+3,a+3) + 1.d0/Hk
              tres(a+3,a+3) = tres(a+3,a+3) + Hkr
              tres(a  ,7) =  nn(a)
              tres(a+3,7) = -nn(a)
              tres(7,a  ) =  nn(a)
              tres(7,a+3) = -nn(a)
            end do ! a

            tres(7,7) = -Ypr

            call invert(tres,7,7)

            do a = 1, 7
              dsol(a) = - tres(a,1)*res(1)
     &                  - tres(a,2)*res(2)
     &                  - tres(a,3)*res(3)
     &                  - tres(a,4)*res(4)
     &                  - tres(a,5)*res(5)
     &                  - tres(a,6)*res(6)
     &                  - tres(a,7)*res(7)
            end do ! a

c         Isotropic hardening only

          else

            do a = 1,3
              tres(a,4) = nn(a)
              tres(4,a) = nn(a)
            end do ! a

            tres(4,4) = -Ypr

            call invert(tres,4,7)

            do a = 1, 4
              dsol(a) = - tres(a,1)*res(1)
     &                  - tres(a,2)*res(2)
     &                  - tres(a,3)*res(3)
     &                  - tres(a,4)*res(7)
            end do ! a
            dsol(7) = dsol(4)
            dsol(4) = 0.0d0
            dsol(5) = 0.0d0
            dsol(6) = 0.0d0
          endif

c         Update Kirchhoff stress and plastic flow

          tau(1) = tau(1) + dsol(1)
          tau(2) = tau(2) + dsol(2)
          tau(3) = tau(3) + dsol(3)
          gam    = gam    + dsol(7)

c         Update accumulated plastic strain

          epp(1) = eppn + sqt23*gam

c         Update Back Stress

          alp(1) = alp(1) + dsol(4)
          alp(2) = alp(2) + dsol(5)
          alp(3) = alp(3) + dsol(6)

c         Update vol.-dev. Kirchhoff stress and stress invariants

          pp     = ( tau(1) + tau(2) + tau(3) ) * one3
          tt(1)  = tau(1) - pp
          tt(2)  = tau(2) - pp
          tt(3)  = tau(3) - pp

          aatr = (alp(1) + alp(2) + alp(3))*one3
          do a = 1,3
            ta(a)   = tt(a) - alp(a) + aatr
          end do ! a

          I1 =  3.0d0 * (pp - aatr)
          J2 = ( ta(1)**2 + ta(2)**2 + ta(3)**2 ) * one2
          J3 = ( ta(1)**3 + ta(2)**3 + ta(3)**3 ) * one3

c         Update vol.-dev. logarithmic strain

          vol_e  = K3inv * pp

          do a = 1, 3
            ee_e(a)  = G2inv * tt(a)
            eps_e(a) = ee_e(a) + vol_e
          end do ! a

          it = it + 1

        end do ! while

c       Warning: check convergence

        if(.not.conv .and. niter.gt.0) then
          write(  *,*) ' *WARNING* No convergence in PLASFD',err,tolc
          write(iow,*) ' *WARNING* No convergence in PLASFD',err,tolc
c         call plstop()
        endif

c       Update elastic left Cauchy-Green tensor and plastic acc. strain

        do a = 1, 3
          if(abs(eps_e(a)).gt.0.0005d0) then
            ll2(a) = exp( 2.0d0 * eps_e(a) ) - 1.0d0
          else
            xe     = 2.d0 * eps_e(a)
            ll2(a) = xe + fac2r*xe**2
     &                  + fac3r*xe**3
     &                  + fac4r*xe**4
     &                  + fac5r*xe**5
     &                  + fac6r*xe**6
          endif
        end do ! a
        do a = 1, 6
          be(a) = ll2(1) * nn_t(p1(a),1)*nn_t(p2(a),1)
     &          + ll2(2) * nn_t(p1(a),2)*nn_t(p2(a),2)
     &          + ll2(3) * nn_t(p1(a),3)*nn_t(p2(a),3)
        end do ! a

c       Update plastic strains

        epl(1) = epl(1) + gam * nn(1)
        epl(2) = epl(2) + gam * nn(2)
        epl(3) = epl(3) + gam * nn(3)

c       Compute elasto-plastic tangent

        do b = 1, 3
          do a = 1, 3
            dtde(a,b) = tres(a,b)
          end do ! a
        end do ! b

      else

c     ELASTIC STEP  ( only tangent computation )

        state  = .false.                   ! Indicate elastic on return

        dtde(1,1) = kk  + four3 * gg
        dtde(1,2) = kk  - two3  * gg
        dtde(1,3) = dtde(1,2)
        dtde(2,1) = dtde(1,2)
        dtde(2,2) = dtde(1,1)
        dtde(2,3) = dtde(1,2)
        dtde(3,1) = dtde(1,2)
        dtde(3,2) = dtde(1,2)
        dtde(3,3) = dtde(1,1)

      end if

c     Save norm of deviator stress term

      epp(2) = sqrt( 2.0d0 * J2 )

c     COMPUTE CAUCHY STRESS

      detr = 1.d0 / detf
      do a = 1, ntm
        sig(a) = (tau(1) * nn_t(p1(a),1)*nn_t(p2(a),1)
     &         +  tau(2) * nn_t(p1(a),2)*nn_t(p2(a),2)
     &         +  tau(3) * nn_t(p1(a),3)*nn_t(p2(a),3)) * detr
      end do ! a
      sig( 9) = epp(1)                   ! Plastic strain for output
      sig(10) = sqrt(1.5d0)*(YY + yield) ! Yield stress value

c     TANGENT TRANSFORMATION

c     Material tangent (computation in the principal basis)

      do a = 1, 3

C       Upper 3x3 block of dd_l()

        do b = 1, 3
          dd_l(b,a) = dtde(b,a)
        end do ! b
        dd_l(a,a) = dd_l(a,a) - 2.0d0 * tau(a)

C       Lower 3x3 block of dd_l() [ diagonal block ]

        b = mod(a,3) + 1
        if (abs(ll2_tr(a)-ll2_tr(b)).gt.tolb) then
          dd_l(a+3,a+3) = ( tau(b) - tau(a)
     &                  +   ll2_tr(a)*tau(b) - ll2_tr(b)*tau(a) ) /
     &                    ( ll2_tr(b) - ll2_tr(a))
        else
          dd_l(a+3,a+3) =  0.5d0*(dtde(a,a) - dtde(b,a)) - tau(a)
        endif

      end do ! a

c     Transform matrix to standard basis

      call tranr4(nn_t,nn_t,ff,.false.)
      call pushr4(ff,ff,dd_l,dd,detf)

      end

      function yfunct(d,epp,I1,J2,J3,f1,f2,f3,
     &                f11,f22,f33,f12,f13,f23,Ypr,YY,gam,ssn,flg)
c    &                f11,f22,f33,f12,f13,f23,Ypr,YY,flg)

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c     2. Add 'gam' to argument list                         06/05/2010
c     3. Add 'ssn' to argument list                         18/05/2010
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'iofile.h'
      include 'pconstant.h'

      real*8  yfunct, yield, gam, ssn

      real*8  d(*),epp(2)
      real*8  I1, J2, J3
      real*8  f1, f2, f3, f11, f22, f33, f12, f13, f23, Ypr
      real*8  J2_12, J2_1, J2_2
      real*8  YY, aa, bb, dd, ss, c1, c2, a1, a2, a3, a4, g1, g2
      real*8  beta, delta, hkin
      logical flg

      integer iyield

      f1  = 0.0d0
      f2  = 0.0d0
      f3  = 0.0d0
      f11 = 0.0d0
      f22 = 0.0d0
      f33 = 0.0d0
      f12 = 0.0d0
      f13 = 0.0d0
      f23 = 0.0d0

      iyield = nint(d(46))           ! Type of yield function
      YY     = d(41) + d(44)*epp(1)  ! sig_y(t_n) yld stress*sqrt(2/3)
      Ypr    = d(44) * sqt23         ! Isotropic yield = 2/3 H_iso
      yield  =-YY                    ! Default return = no yield

c     VON MISES        YIELD FUNCTION

      if (iyield.eq.1) then

        aa    = (d(41) - d(42))*exp(-d(43)*epp(1))
        YY    = d(42)  + aa + d(44)*epp(1)
        Ypr   = Ypr    - sqt23*d(43)*aa
        ss    = sqrt( 2.0d0 * J2 )
        yield = ss - YY

c       Compute yield function derivatives

        if (flg) then

          if (ss.ne.0.0d0) then
            f2  =   1.0d0 / ss
            f22 = - f2**3
          end if

        end if

c     DRUCKER PRAGER   YIELD FUNCTION

      elseif (iyield.eq.2) then

        aa    = d(42)
        ss    = sqrt( 2.0d0 * J2 )
        yield = ss + one3 * aa * I1 - YY

c       Compute yield function derivatives

        if (flg) then

          f1  = one3 * aa

          if (ss.ne.0.0d0) then
            f2  =   1.0d0 / ss
            f22 = - f2**3
          end if

        end if

c     PRAGER-LODE      YIELD FUNCTION

      elseif (iyield.eq.3) then

        if (J2.gt.0.0d0) then

          if (flg) then

c           Read in material parameters

            bb      = d(42)

c           Compute yield function

            c1      = 1.0d0 / sqrt(2.0d0)
            c2      = sqrt(13.5d0) * bb

c           yield   = J2 * ( 1.0d0 + bb * J3 * J2 **(-1.5d0)) - YY
            yield   = ( sqrt(2.0d0*J2) + c2 * J3 / J2 ) - YY

c           Compute yield function derivatives

            J2_1    =   1.d0/J2
            J2_2    =   J2_1**2
            J2_12   =   sqrt(J2_1)

            f2      =   c1 * J2_12            - c2 * J3 * J2_2
            f3      =                           c2      * J2_1
            f22     = - c1 * J2_12**3 * 0.5d0 + c2 * J3 * J2_1**3 * 2.d0
            f23     =                         - c2      * J2_2
          end if

        else
          yield = -YY
        end if

c     GENERALIZED MISES PLASTICITY

      elseif (iyield.eq.4) then

        ss = sqrt( 2.0d0 * J2 )           ! |S_tr|
        a1 = ss - YY                      ! |S_tr|-sqrt(2/3)*sig_y(t_n)
        a2 = ss - ssn                     ! |S_tr| - |S_n|
        if(min(a1,a2).gt.0.0d0) then
          beta  =  d(42) - d(41)          ! beta
          delta =  d(43)                  ! delta * 2/3
          hkin  =  d(45)                  ! H_kin * 2/3
          if(.not.flg) then
            a3    =  2.d0*d(22) - delta       ! 2*G - delta
            a4    =  delta      + Ypr + hkin  ! delta+2/3(H_iso+H_kin)
            g1    =  2.d0*d(22) + hkin        ! 2*G  +2/3 H_kin
            g2    =  g1 + Ypr                 ! 2*G  +2/3(H_iso+H_kin)

c           Coefficients of quadratic equation

            aa = g2*a3
            bb = a1*a3 + a2*g2 + a4*beta

c           Compute discriminant

            dd = bb*bb - 4.d0*aa*a1*a2
            if(dd.lt.0.0d0) then
              write(iow,*) ' Gen Plastic Discriminant Error =',dd
            endif
            gam   = 0.5d0 *(bb - sqrt(abs(dd)))/aa

c           Update functions

            ss = ss - g1*gam                     ! |SS_n+1|
            YY = YY + d(44)*gam*sqt23            ! sqrt(2/3)*sig_y,n+1
          endif

          a1 = ss - YY
          a2 = ss - ssn + hkin*gam            ! epp(2) = |SS_n|
          bb = (a1*hkin - (a2 + delta*gam)*Ypr
     &       -  delta*(beta - a1) - beta*(hkin + Ypr))            ! RLT
     &       / (a1 + a2 + delta*gam)
c    &       -  delta*(beta - a1 - Ypr*gam) - beta*(hkin + Ypr))  ! FDA

c         Return Values

          Ypr = -bb                     ! Makes diagonal for lambda = bb
          if(flg) then
            yield = 0.0d0
            if(ss.ne.0.0d0) then
              f2  =   1.0d0/ss
              f22 = - f2**3
            endif
          else
            yield = (YY + Ypr*gam)/(a1 + a2 + delta*gam)
          endif
        else
          gam   = 0.0d0
        endif

c     NO YIELD FUNCTION SPECIFIED

      else

        yield = 0.0d0
        write (iow,9000)
        if (ior.lt.0) then
          write (*,9000)
        endif

      end if

      yfunct = yield

c     Format

9000  format(' *** NO YIELD FUNCTION SPECIFIED ***')

      end
