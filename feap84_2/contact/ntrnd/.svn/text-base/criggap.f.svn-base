c$Id:$
      subroutine criggap(cs0,ch2,cn,xs,ndm, rsd,tan, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Extract transformation data for cs0 to position  15/05/2013
c          contact rigid cylindrical and spherical surface
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute active elements and residual & tangent for node
c              to rigid surface contacts.

c              Gap: gn = nn*(sqrt(xs(i)*xs(i)) - rr - un*prop(np))

c     Inputs:
c       cs0(*,*)  - Surface properties:
c              nsopt = cs0(2,0): 1 = cylinder, 2 = sphere, 3 = plane
c              dn    = cs0(3,0): Normal direction from rigid surface
c                                     (1 <= dn <= ndm)
c              rr    = Radius of cylinder or sphere; Location of plane
c              un    = Displacement amplitude for surface movement
c              np    = Proportional load number
c       ch2(*)    - History variables
c       cn(*)     - Penalty parameters: k0 = cn(1): Linear constant
c                                       k3 = cn(2): Cubic  constant
c       xs(3,*)   - Coordinates for contact node (slave)
c       isw       - 1: Compute active element;
c                   2: Compute residual & tangent

c     Outputs:
c       rsd(4)    - Residual array
c       tanm(4,4) - Tangent array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_comnd.h'
      include   'c_keyh.h'
      include   'c_pair.h'
      include   'counts.h'
      include   'prld1.h'
      include   'prlod.h'

      integer    ndm, isw, i,j, np, nsopt, dn
      real*8     cs0(nr0,n0c1:*),ch2(*),xs(3),rsd(4),tan(4,4)
      real*8     xx(3), r(3), rl(3), dtn(4)
      real*8     gn, rr,un, fn, cn(*), xl, ri, fl, nn

      real*8     tr(3,3),xr(3)

c     Extract surface type and parameters

      nsopt = nint(cs0(2,0))
      dn    = nint(cs0(3,0))
      rr    = cs0(4,0)
      un    = cs0(5,0)
      np    = nint(cs0(6,0))

c     Extract transformation array

      tr(1,1) = cs0( 7,0)
      tr(2,1) = cs0( 8,0)
      tr(3,1) = cs0( 9,0)
      tr(1,2) = cs0(10,0)
      tr(2,2) = cs0(11,0)
      tr(3,2) = cs0(12,0)
      tr(1,3) = cs0(13,0)
      tr(2,3) = cs0(14,0)
      tr(3,3) = cs0(15,0)
      xr(1)   = cs0(16,0)
      xr(2)   = cs0(17,0)
      xr(3)   = cs0(18,0)

c     Set displacement value of surface

      if(np.gt.0) then
        un = un*prldv(np)
      else
        un = un*prop
      endif

c     Zero residual and tangent

      if(isw.eq.2) then
        do i = 1,4
          rsd(i) = 0.0d0
          dtn(i) = 1.0d0
          do j = 1,4
            tan(j,i) = 0.0d0
          end do ! j
        end do ! i
      endif

c     Cylindrical or Spherical surface

      if(nsopt.le.2) then

        do i = 1,3
          xx(i) = tr(1,i)*(xs(1) - xr(1))
     &          + tr(2,i)*(xs(2) - xr(2))
     &          + tr(3,i)*(xs(3) - xr(3))
        end do ! i
        if(nsopt.eq.1) then
          xx(3) = 0.0d0
          do i = 3,4
            dtn(i) = 0.0d0
          end do ! i
        endif

c       Compute gap function

        xl = sqrt(xx(1)**2 + xx(2)**2 + xx(3)**2)
        if(dn.gt.0) then
          gn =  xl - rr - un
          nn =  1.d0
        else
          gn =  rr + un - xl
          nn = -1.d0
        endif
        ch2(p1(2)) = gn
        if(ifaugm.le.1) then
          fn = (cn(1) + cn(2)*gn*gn)*gn
        else
          fn = (cn(1) + cn(2)*gn*gn)*gn + ch2(p1(151))
        endif

c       Lagrange multiplier

        if(ifsolm.eq.2) then
          fn = fn + ch2(p1(21))
        endif

c       Test on contact force

        if(isw.eq.1) then

c         Set contact state

          if(fn.lt.0.0d0) then
            ch2(p1(1)) = 1.d0
          elseif(gn.lt.0.0d0) then
            ch2(p1(1)) = 1.d0
          else
            ch2(p1(1)) = 0.d0
          endif

c       Compute residual and tangent arrays

        else

c         Compute local normal

          do i = 1,3
            rl(i) = nn*xx(i)/xl
          end do ! i

c         Transform to global frame

          r(:) = (tr(:,1)*rl(1) + tr(:,2)*rl(2) + tr(:,3)*rl(3))

c         Compute global tangent matrix

          if(niter.gt.0) then
            fl = fn/xl*nn
          else
            fl = 0.0d0
          endif
          do i = 1,ndm
            ri = (cn(1) + 3.d0*gn*gn*cn(2) - fl)*r(i)
            do j = 1,ndm
              tan(i,j) = ri*r(j)
            end do ! j
            tan(i,i) = tan(i,i) + fl*dtn(i)
          end do ! i

c         Compute global residual matrix

          do i = 1,ndm
            rsd(i) = -r(i)*fn
          end do ! i

c         Lagrange multiplier terms

          if(ifsolm.eq.2) then
            rsd(4) = -gn
            do i = 1,ndm
              tan(i,4) = r(i)
              tan(4,i) = r(i)
            end do ! i
            tan(4,4) = 0.0d0
          endif

        endif ! isw

c     Cartesian plane surface

      elseif(nsopt.eq.3) then

        if(dn.gt.0) then
          gn  = xs(dn) - rr - un
          nn  = 1.d0
        else
          gn = rr + un - xs(-dn)
          nn = -1.d0
        endif
        ch2(p1(2)) = gn
        if(ifaugm.le.1) then
          fn = (cn(1) + cn(2)*gn*gn)*gn
        else
          fn = (cn(1) + cn(2)*gn*gn)*gn + ch2(p1(151))
        endif

c       Lagrange multiplier

        if(ifsolm.eq.2) then
          fn = fn + ch2(p1(21))
        endif

c       Test on contact force

        if(isw.eq.1) then

c         Set contact state

          if(fn.lt.0.0d0) then
            ch2(p1(1)) = 1.d0
          else
            ch2(p1(1)) = 0.d0
          endif

c       Compute residual and tangent arrays

        else

c         Compute global residual and tangent

          dn         =  abs(dn)
          rsd(dn)    = -fn*nn
          tan(dn,dn) =  cn(1) + 3.d0*gn*gn*cn(2)

c         Lagrange multiplier terms

          if(ifsolm.eq.2) then
            rsd(4)    = -gn
            tan(dn,4) =  nn
            tan(4,dn) =  nn
          endif

        endif ! isw

      endif ! nsopt

      end
