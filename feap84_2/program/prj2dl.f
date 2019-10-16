c$Id:$
      subroutine prj2dl(gap0,tolg,tol,x0,ip,xin,x,xl,ndm,numnp,polfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded by:                           Date:                   rel.:
c         Robert L. Taylor                 November 1, 1998          1.0

c      Purpose: Project points onto line segment.  Line segment may
c               be represented in polar or cartesian coordinates.

c      Inputs:
c         gap0(2) - Specified gap tolerances
c         tolg    - Specified gap tolerance
c         tol     - Specified search tolerance
c         x0(*)   - Center of polar frame
c         x(*)    - Nodal coordinates
c         xl(*)   - Patch polar coordinates
c         ndm     - Space dimension of mesh
c         numnp   - Number of nodes in mesh
c         polfl   - Polar flag

c      Scratch:
c         ip(*)   - Nodal integer list storage
c         xin(*)  - Nodal real list storage

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pconstant.h'

      logical   polfl
      integer   i,n,nits, ndm,numnp
      integer   ip(numnp)
      real*8    tolg,tol,tolxi,d,gap,gap0(2),th
      real*8    dx,dy,dz,ff,df,xi,dxi, xin(*),x(ndm,*),xmin(3),xmax(3)
      real*8    xx(3),xp(3),x0(3),xl(3,3),dxl(3),dxq(3),vv(3)

      save

      data      tolxi / 1.0d-8 /
      data      nits /50/

c     Set angle degree to radian

      th     = 180.d0/pi

c     Find maximum and minimum of search region

      do i = 1,3
        xmin(i) = min(xl(i,1),xl(i,2),xl(i,3))
        xmax(i) = max(xl(i,1),xl(i,2),xl(i,3))
      end do ! i

c     Set gap condition

      if(gap0(1).gt.0.0d0) then
        gap = gap0(1)
      else
        gap = 0.5d0*tolg*max(xmax(1)-xmin(1),xmax(2)-xmin(2),
     &                       xmax(3)-xmin(3))
      endif

      if(polfl) then
        ff      = (xmax(1) - xmin(1))*gap + 1.d-10*(xmax(1) + xmin(1))
     &          + gap0(2)
        xmin(1) = xmin(1) - 1.5d0*ff
        xmax(1) = xmax(1) + 1.5d0*ff
        xmax(2) = xmax(2) + gap
        xmin(2) = xmin(2) - gap
        vv(1)   = cos(xl(2,3)/th)
        vv(2)   = sin(xl(2,3)/th)
        ff      = (xmax(3) - xmin(3))*gap
        xmin(3) = xmin(3) - 1.5d0*ff
        xmax(3) = xmax(3) + 1.5d0*ff
      else
        th = max(xmax(1)-xmin(1),xmax(2)-xmin(2),xmax(3)-xmin(3))
        do i = 1,3
          xmin(i) = xmin(i) - (gap0(2) + 0.1d0)*th
          xmax(i) = xmax(i) + (gap0(2) + 0.1d0)*th
        end do ! i
      endif

c     Tag nodes within surface loading area

      do n = 1,numnp
        ip(n)  = 0
        xin(n) = 0.0d0
      end do ! n

c     Find nodes on surface (within gap)

      do i = 1,3
        dxl(i) = 0.5d0*(xl(i,2) - xl(i,1))
        dxq(i) = xl(i,1) + xl(i,2) - 2.d0*xl(i,3)
      end do ! i

      do n = 1,numnp
        if(polfl) then
          dx    =  vv(1)*(x(1,n) - x0(1)) + vv(2)*(x(2,n) - x0(2))
          dy    = -vv(2)*(x(1,n) - x0(1)) + vv(1)*(x(2,n) - x0(2))
          xx(2) = atan2(dy,dx)*th + xl(2,3)
          xx(1) = sqrt((x(1,n) - x0(1))**2 + (x(2,n) - x0(2))**2)
          if(ndm.eq.3) then
            xx(3) = x(3,n) - x0(3)
          else
            xx(3) = 0.0d0
          endif
        else
          xx(1) = x(1,n)
          xx(2) = x(2,n)
          if(ndm.eq.3) then
            xx(3) = x(3,n)
          else
            xx(3) = 0.0d0
          endif
        endif
        if((xx(1).ge.xmin(1) .and. xx(1).le.xmax(1)) .and.
     &     (xx(2).ge.xmin(2) .and. xx(2).le.xmax(2)) .and.
     &     (xx(3).ge.xmin(3) .and. xx(3).le.xmax(3)) ) then

c         Compute location of closest point on surface

          xi = -1.0d0
          ff = (xl(1,1) - xx(1))**2
     &       + (xl(2,1) - xx(2))**2
     &       + (xl(3,1) - xx(3))**2
          do i = 2,3
            df = (xl(1,i) - xx(1))**2
     &         + (xl(2,i) - xx(2))**2
     &         + (xl(3,i) - xx(3))**2
            if(df.lt.ff) then
              ff = df
              xi = dble(3-i)
            endif
          end do ! i
          do i = 1,nits
            xp(1) = xl(1,3) + xi*(dxl(1) + 0.5d0*xi*dxq(1))
            xp(2) = xl(2,3) + xi*(dxl(2) + 0.5d0*xi*dxq(2))
            xp(3) = xl(3,3) + xi*(dxl(3) + 0.5d0*xi*dxq(3))
            dx    = dxl(1)  + xi*dxq(1)
            dy    = dxl(2)  + xi*dxq(2)
            dz    = dxl(3)  + xi*dxq(3)
            ff    = dx*(xp(1) - xx(1))
     &            + dy*(xp(2) - xx(2))
     &            + dz*(xp(3) - xx(3))
            df    = dx*dx + dy*dy + dz*dz
            if(i.gt.3) then
              df = df + dxq(1)*(xp(1) - xx(1))
     &                + dxq(2)*(xp(2) - xx(2))
     &                + dxq(3)*(xp(3) - xx(3))
            endif
            dxi   = -ff/df
            xi    = xi + dxi
            if(abs(dxi).lt.tolxi) go to 200
          end do ! i

c         Non linear interations did not converge

          write(  *,*) ' No convergence in PRJ2DL. Node =',n,dxi
          write(iow,*) ' No convergence in PRJ2DL. Node =',n,dxi

c         Compute gap and if too large eliminate this node

200       d = sqrt((xx(1) - xp(1))**2
     &           + (xx(2) - xp(2))**2
     &           + (xx(3) - xp(3))**2)
c         if(d .lt. gap .and. abs(xi).lt.3.0d0+tol) then
          if(d .lt. gap .and. abs(xi).lt.100.0d0+tol) then
            xin(n) = xi
            ip(n)  = 1
          else
            ip(n) = 0
          endif
        endif
      end do ! n

      end
