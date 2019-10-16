c$Id:$
      subroutine modlf_apl(d, detf, f, epn, istrt, ds, sig, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    12/12/2012
c       1 Add 'istrt' option                                20/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Orthotropic elastic-plasticity in logarithmic strain
c               strain energy function

c      Inputs:
c        d(*)     - Parameter array
c        detf     - Determinant of deformation gradient at t_n+1
c        f(3,3)   - Deformation gradient at t_n+1
c        epn(*)   - Current plastic strains
c        istrt    - Plasticity control for first iteration
c        isw      - Switch parameter: isw=14 initialize arrays.

c      Outputs:
c        ds(6,6)  - Tangent moduli
c        sig(6)   - Cauchy stress
c-----[--.----+----.----+----.-----------------------------------------]
c     w = 1/2 k(tre)^2 + 1/2 mu_1 (e_11-e_22)^2 + 3/2 mu_2 (e'_33)^2
c       + 2mu_3(e_12)^2 + 2mu_4(e_23)^2 + 2mu_5(e_13)^2
c       + beta_1(e_11-e_22)(e'_33)+beta_5(e_11-e_22)tre +beta_6e'_33tre
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pconstant.h'

      logical    plast,orth
      integer    istrt, isw, i,j
      real*8     detf, temp
      real*8     d(*),f(3,3),rf(3,3),ds(6,6),sig(13),epn(*)
      real*8     ss(6),ro(3,3),br(3),ee(6), dp(9), sp(6)
      real*8     t(6,6),gs(6,6),dd(6,6)
      real*8     a2(3), a3(3)
      real*8     sr(3),fr(6),dr(3),ddr(3),dti(6)

      real*8     rr(3,3)
      integer    rot

      save

      data       a2  /0.0d0, 1.0d0, -1.0d0/

c     Initialize history variables

      if(isw.eq.14) then
        epn(1:14) = 0.d0
        return
      endif

c     Set a3

      a3(1) =  two3
      a3(2) = -one3
      a3(3) = -one3

c     Rotate coordinate to align with symmetry axises

      orth = nint(d(242)).eq.1
      if(orth) then
        call vectf_apl(d(243),f, rf)
      else
        rf(:,:) = f(:,:)
      endif

c     Calculate Green deformation tensor.

      ro(1,1) = rf(1,1)*rf(1,1) + rf(2,1)*rf(2,1) + rf(3,1)*rf(3,1)
      ro(2,2) = rf(1,2)*rf(1,2) + rf(2,2)*rf(2,2) + rf(3,2)*rf(3,2)
      ro(3,3) = rf(1,3)*rf(1,3) + rf(2,3)*rf(2,3) + rf(3,3)*rf(3,3)

      ro(1,2) = rf(1,1)*rf(1,2) + rf(2,1)*rf(2,2) + rf(3,1)*rf(3,2)
      ro(2,1) = ro(1,2)

      ro(2,3) = rf(1,2)*rf(1,3) + rf(2,2)*rf(2,3) + rf(3,2)*rf(3,3)
      ro(3,2) = ro(2,3)

      ro(3,1) = rf(1,3)*rf(1,1) + rf(2,3)*rf(2,1) + rf(3,3)*rf(3,1)
      ro(1,3) = ro(3,1)

      call eig3(ro,br,rot)

c     Set logarithmic strain measure relations

      sr(1:3)  = 0.5d0                       ! f'*f
      fr(1:3)  = 0.5d0*log(br(1:3))          ! f
      dr(1:3)  = 0.5d0/br(1:3)               ! f'
      ddr(1:3) =-dr(1:3)/br(1:3)             ! f''
      fr(4:6)  = 0.0d0

c     Eulerian triad

      do j = 1,3
        temp = 1.0d0/sqrt(br(j))
        rr(:,j) = (rf(:,1)*ro(1,j)
     &          +  rf(:,2)*ro(2,j)
     &          +  rf(:,3)*ro(3,j))*temp
      end do ! j

      call pushr2(ro,fr,ee,1.0d0)

c     Adjust for any existing plastic strain

      plast  =  nint(d(40)).eq.6
      if(plast) then
        ee(1:6) = ee(1:6) - epn(1:6)
      end if

c     Stress response

      ds = 0.0d0

c     Compute Papdopoulos-Lu elastic parameters from input parameters

      dp(1) = (d(21) + d(22) + d(23)               ! K
     &      + (d(24) + d(25) + d(26))*2.d0)*one9
      dp(2) = (d(22) + d(23) - 2.d0*d(25))*0.25d0  ! mu_2
      dp(3) = (d(21) - d(24) - d(26))*one3         ! mu_3
     &      + (d(22) + d(23) + 2.0d0*d(25))*one12
      dp(4) = (d(24) - d(26))*0.25d0               ! beta_1
     &      + (d(23) - d(22))*0.125d0
      dp(5) = (d(24) - d(26) + d(22) - d(23))*one6 ! beta_5
      dp(6) = (d(21) - d(25))*one3                 ! beta_6
     &      + (d(24) + d(26) - d(22) - d(23))*one6
      dp(7) =  d(27)                               ! mu_4
      dp(8) =  d(28)                               ! mu_5
      dp(9) =  d(29)                               ! mu_6

c     Elastic response

      do j = 1, 3
        ds(1:3,j) = dp(1) + dp(2)*a2(:)*a2(j)
     &            + 3.0d0 * dp(3)*a3(:)*a3(j)
     &            + dp(4) * (a2(:)*a3(j) + a3(:)*a2(j))
     &            + dp(5) * (a2(:)+a2(j)) + dp(6)*(a3(:)+a3(j))
      end do ! j
      ds(4,4) = dp(7)
      ds(5,5) = dp(8)
      ds(6,6) = dp(9)

c     Stress response: Can generalize for rotated elastic moduli

        ss(:) = ds(:,1)*ee(1) + ds(:,2)*ee(2) + ds(:,3)*ee(3)
     &        +(ds(:,4)*ee(4) + ds(:,5)*ee(5) + ds(:,6)*ee(6))*2.d0

      if(plast) then
        call plasf_apl(d,dp,epn,istrt,ss,ee,ds)
        sig(7) = epn(7)
      end if

c     Transpose 'ro'

      do j = 1,2
        do i = j+1,3
          temp = ro(i,j)
          ro(i,j) = ro(j,i)
          ro(j,i) = temp
        end do ! i
      end do ! j

c     Stress component in principal direction

      call pushr2(ro,ss,sp,1.0d0)

c     Geometric stiffness

      call geomf_apl(t,rr,ro,sp,gs, br,sr,fr,dr,ddr)

c     Transform Cauchy stress to standard basis

      temp     = 1.0d0/detf
      sig(1:6) = (t(:,1)*ss(1) + t(:,2)*ss(2) + t(:,3)*ss(3)
     &         +  t(:,4)*ss(4) + t(:,5)*ss(5) + t(:,6)*ss(6))*temp

c     Tranform moduli to standard basis

      dd = ds
      do i = 1,6

        dti(:) = t(i,1)*dd(1,:) + t(i,2)*dd(2,:) + t(i,3)*dd(3,:)
     &         + t(i,4)*dd(4,:) + t(i,5)*dd(5,:) + t(i,6)*dd(6,:)

c       Compute spatial tensor: ds = dti*t_tran

        ds(i,:) = (dti(1)*t(:,1) + dti(2)*t(:,2) + dti(3)*t(:,3)
     &          +  dti(4)*t(:,4) + dti(5)*t(:,5) + dti(6)*t(:,6)
     &          +  gs(i,:))*temp

      end do ! i

      end
