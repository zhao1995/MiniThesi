c$Id:$
      subroutine geomf_apl(tt,rr,ro,ss,ds, br,sr,fr,dr,ddr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    12/12/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Geometric stiffness for finite orthotropic plasticity.
c               Calculate 1st & 2nd strain rate transformation matrix
c               in principal direction.

c      Inputs:
c        rr(3,3)   - Eulerian triad
c        ro(3,3)   - Green deformation tensor
c        ss(6)     - Green stress
c        ds(6,6)   - Material tangent moduli
c        br(3)     - lambda^2
c        sr(3)     - F' * F
c        fr(3)     - F
c        dr(3)     - F'

c      Outputs:
c        tt(6,6)   - Transformation array
c        ds(6,6)   - Material + Geometric tangent moduli
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pconstant.h'

      integer    i
      real*8     gs(6,6),ss(6),rr(3,3),ro(3,3),t(6,6),f(6,6)
      real*8     ds(6,6),ff(6,6),tt(6,6),dti(6),dsi(6)
      real*8     br(3),sr(3),fr(6),dr(3),ddr(3)
      real*8     fac1,fac2,fac3,rt12,rt23,rt31,tol

      data       tol/1.0d-6/

      gs = 0.0d0
      ff = 0.0d0

c     Default values at three repeated eigenvalues

      do i = 1, 3
        gs(i,i) = 4.0d0*ss(i)*ddr(i)*br(i)*br(i)
        ff(i,i) = 2.0d0*sr(i)
      end do ! i
      gs(4:6,4:6) = 0.5d0*gs(1:3,1:3)
      ff(4:6,4:6) = ff(1:3,1:3)

      gs(1,4) = 2.0d0*ss(4)*ddr(1)*br(1)*br(1)
      gs(2,4) = 2.0d0*ss(4)*ddr(1)*br(1)*br(1)
      gs(2,5) = 2.0d0*ss(5)*ddr(2)*br(2)*br(2)
      gs(3,5) = 2.0d0*ss(5)*ddr(2)*br(2)*br(2)
      gs(1,6) = 2.0d0*ss(6)*ddr(3)*br(3)*br(3)
      gs(3,6) = 2.0d0*ss(6)*ddr(3)*br(3)*br(3)
      gs(5,6) =       ss(4)*ddr(1)*br(1)*br(1)
      gs(4,6) =       ss(5)*ddr(2)*br(2)*br(2)
      gs(4,5) =       ss(6)*ddr(3)*br(3)*br(3)

      rt12 = sqrt(br(1)*br(2))
      rt23 = sqrt(br(2)*br(3))
      rt31 = sqrt(br(3)*br(1))

c     Modification in case of repeated eigenvalues

      if( abs(br(1) - br(2)).gt.tol) then
        ff(4,4) = 2.0d0*(fr(1)-fr(2))/(br(1)-br(2))*rt12
        fac1    = (fr(2)-fr(1)-dr(1)*(br(2)-br(1)))/(br(2)-br(1))**2
        fac2    = (fr(1)-fr(2)-dr(2)*(br(1)-br(2)))/(br(1)-br(2))**2
        gs(4,4) = 2.0d0*(ss(1)*fac1+ss(2)*fac2)*br(1)*br(2)
        gs(1,4) = 4.0d0*ss(4)*fac1*br(1)*rt12
        gs(2,4) = 4.0d0*ss(4)*fac2*br(2)*rt12
      end if

      if( abs(br(2) - br(3)).gt.tol) then
        ff(5,5) = 2.0d0*(fr(2)-fr(3))/(br(2)-br(3))*rt23
        fac1    = (fr(3)-fr(2)-dr(2)*(br(3)-br(2)))/(br(3)-br(2))**2
        fac2    = (fr(2)-fr(3)-dr(3)*(br(2)-br(3)))/(br(2)-br(3))**2
        gs(5,5) = 2.0d0*(ss(2)*fac1+ss(3)*fac2)*br(2)*br(3)
        gs(2,5) = 4.0d0*ss(5)*fac1*br(2)*rt23
        gs(3,5) = 4.0d0*ss(5)*fac2*br(3)*rt23
      end if

      if( abs(br(3) - br(1)).gt.tol) then
        ff(6,6) = 2.0d0*(fr(3)-fr(1))/(br(3)-br(1))*rt31
        fac1    = (fr(1)-fr(3)-dr(3)*(br(1)-br(3)))/(br(1)-br(3))**2
        fac2    = (fr(3)-fr(1)-dr(1)*(br(3)-br(1)))/(br(3)-br(1))**2
        gs(6,6) = 2.0d0*(ss(3)*fac1+ss(1)*fac2)*br(3)*br(1)
        gs(3,6) = 4.0d0*ss(6)*fac1*br(3)*rt31
        gs(1,6) = 4.0d0*ss(6)*fac2*br(1)*rt31
      end if

c     LL_ijjkkl terms

      if    ( abs(br(1)-br(2)).gt.tol.and.abs(br(2)-br(3)).gt.tol
     &   .and.abs(br(3)-br(1)).gt.tol) then
        fac1 = (fr(1)+fr(2)-2.0d0*fr(3))/(br(2)-br(3))/(br(3)-br(1))
     &       + (fr(1)-fr(2))/(br(1)-br(2))/(br(2)-br(3))
     &       - (fr(1)-fr(2))/(br(1)-br(2))/(br(3)-br(1))
        fac2 = (fr(2)+fr(3)-2.0d0*fr(1))/(br(3)-br(1))/(br(1)-br(2))
     &       + (fr(2)-fr(3))/(br(2)-br(3))/(br(3)-br(1))
     &       - (fr(2)-fr(3))/(br(2)-br(3))/(br(1)-br(2))
        fac3 = (fr(3)+fr(1)-2.0d0*fr(2))/(br(1)-br(2))/(br(2)-br(3))
     &       + (fr(3)-fr(1))/(br(3)-br(1))/(br(1)-br(2))
     &       - (fr(3)-fr(1))/(br(3)-br(1))/(br(2)-br(3))
        gs(5,6) = fac1*ss(4)*br(3)*rt12
        gs(4,6) = fac2*ss(5)*br(1)*rt23
        gs(4,5) = fac3*ss(6)*br(2)*rt31
      else if(abs(br(1)-br(2)).le.tol.and.abs(br(2)-br(3)).gt.tol
     &   .and.abs(br(3)-br(1)).gt.tol) then
        fac1    = 2.d0*(fr(3)-fr(1)-dr(1)
     &                *(br(3)-br(1)))/(br(3)-br(1))**2
        gs(5,6) = fac1*ss(4)*br(3)*rt12
        gs(4,6) = fac1*ss(5)*br(1)*rt23
        gs(4,5) = fac1*ss(6)*br(2)*rt31
      else if(abs(br(1)-br(2)).gt.tol.and.abs(br(2)-br(3)).le.tol
     &   .and.abs(br(3)-br(1)).gt.tol) then
        fac1    = 2.d0*(fr(1)-fr(2)-dr(2)
     &                *(br(1)-br(2)))/(br(1)-br(2))**2
        gs(5,6) = fac1*ss(4)*br(3)*rt12
        gs(4,6) = fac1*ss(5)*br(1)*rt23
        gs(4,5) = fac1*ss(6)*br(2)*rt31
      else if(abs(br(1)-br(2)).gt.tol.and.abs(br(2)-br(3)).gt.tol
     &     .and.abs(br(3)-br(1)).le.tol) then
        fac1    = 2.d0*(fr(2)-fr(3)-dr(3)
     &                *(br(2)-br(3)))/(br(2)-br(3))**2
        gs(5,6) = fac1*ss(4)*br(3)*rt12
        gs(4,6) = fac1*ss(5)*br(1)*rt23
        gs(4,5) = fac1*ss(6)*br(2)*rt31
      end if

c     Symmetric part

      gs(4,1) = gs(1,4)
      gs(4,2) = gs(2,4)
      gs(5,2) = gs(2,5)
      gs(5,3) = gs(3,5)
      gs(6,1) = gs(1,6)
      gs(6,3) = gs(3,6)
      gs(6,5) = gs(5,6)
      gs(6,4) = gs(4,6)
      gs(5,4) = gs(4,5)

c     Push to another coordinate system

      call strnf_apl(rr,t)
      call strnf_apl(ro,f)

c     t*gs*t', t*ff*f

      do i = 1,6
          dti(:) = t(i,1)*ff(1,:) + t(i,2)*ff(2,:) + t(i,3)*ff(3,:)
     &           + t(i,4)*ff(4,:) + t(i,5)*ff(5,:) + t(i,6)*ff(6,:)
          dsi(:) = t(i,1)*gs(1,:) + t(i,2)*gs(2,:) + t(i,3)*gs(3,:)
     &           + t(i,4)*gs(4,:) + t(i,5)*gs(5,:) + t(i,6)*gs(6,:)

          tt(i,:) = dti(1)*f(1,:) + dti(2)*f(2,:) + dti(3)*f(3,:)
     &            + dti(4)*f(4,:) + dti(5)*f(5,:) + dti(6)*f(6,:)
          ds(i,:) = dsi(1)*t(:,1) + dsi(2)*t(:,2) + dsi(3)*t(:,3)
     &            + dsi(4)*t(:,4) + dsi(5)*t(:,5) + dsi(6)*t(:,6)

      end do ! i

      end
