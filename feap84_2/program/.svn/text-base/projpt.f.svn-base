c$Id:$
      integer function projpt(x,xp,xi,gap,normal,shp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct geometric tangent term                   14/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Normal to 9-node patch

c      Inputs:
c         x(ndm,*)       - Nodal coordinates
c         xp(ndm)        - Element coordinates
c         xi             - Natural coordinate
c         gap            - Specified gap

c      Outputs:
c         normal(3)      - Normal vector
c         shp(9)         - Shape function
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   noconv
      integer   i,n,iters
      real*8    gap,detr, tol, a11,a12,a22
      real*8    r(2), t1(3),t2(3), dt11(3),dt12(3),dt22(3)
      real*8    dx(3),xi(2),dxi(2),shp(9),d1shp(2,9),d2shp(3,9)
      real*8    normal(3),x(3),xp(3,9)

      data      tol /1.d-6/

c     Projections to nine node patch

      noconv = .true.
      iters  = 0
      xi(1)  = 0.0d0
      xi(2)  = 0.0d0
      do while(noconv)
        iters = iters + 1
        call pshp9(xi,shp,d1shp,d2shp)
        do i = 1,3
          dx(i)   = -x(i)
          t1(i)   = 0.0d0
          t2(i)   = 0.0d0
          dt11(i) = 0.0d0
          dt12(i) = 0.0d0
          dt22(i) = 0.0d0
          do n = 1,9
            dx(i)   = dx(i)   + shp(n)*xp(i,n)
            t1(i)   = t1(i)   + d1shp(1,n)*xp(i,n)
            t2(i)   = t2(i)   + d1shp(2,n)*xp(i,n)
            dt11(i) = dt11(i) + d2shp(1,n)*xp(i,n)
            dt12(i) = dt12(i) + d2shp(2,n)*xp(i,n)
            dt22(i) = dt22(i) + d2shp(3,n)*xp(i,n)
          end do ! n
        end do ! i

c       Form residual

        r(1) = dx(1)*t1(1) + dx(2)*t1(2) + dx(3)*t1(3)
        r(2) = dx(1)*t2(1) + dx(2)*t2(2) + dx(3)*t2(3)

c       Form linearized tangent terms

        a11  = t1(1)*t1(1) + t1(2)*t1(2) + t1(3)*t1(3)
        a12  = t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3)
        a22  = t2(1)*t2(1) + t2(2)*t2(2) + t2(3)*t2(3)

c       Add geometric stiffness terms

        if(iters.gt.3) then
          a11 = a11
     &        + dx(1)*dt11(1) + dx(2)*dt11(2) + dx(3)*dt11(3)
          a12 = a12
     &        + dx(1)*dt12(1) + dx(2)*dt12(2) + dx(3)*dt12(3)
          a22 = a22
     &        + dx(1)*dt22(1) + dx(2)*dt22(2) + dx(3)*dt22(3)
        endif

c       Solve Newton equations

        detr   = 1.d0/(a11*a22 - a12*a12)
        dxi(1) = (-a22*r(1) + a12*r(2))*detr
        dxi(2) = ( a12*r(1) - a11*r(2))*detr

        xi(1)  = xi(1) + dxi(1)
        xi(2)  = xi(2) + dxi(2)

        if((max(abs(dxi(1)),abs(dxi(2))).lt.tol)
     &      .or. (iters .gt. 100)              ) noconv = .false.

      end do ! while

c     Stop on non-convergence

      if(iters.gt.100) then
        write(ilg,3000) 'No convergence.'
        write(iow,3000) 'No convergence.'
        call plstop()
      endif

c     Check position of projection: projpt = 1 "on" surface; = 1 "not"

      if((dx(1)**2+dx(2)**2+dx(3)**2 .lt. gap*gap) .and.
     &   (max(abs(xi(1)),abs(xi(2))) .lt. 1.d0+1.d2*tol) ) then
        projpt    = 1
        normal(1) = t1(2)*t2(3) - t1(3)*t2(2)
        normal(2) = t1(3)*t2(1) - t1(1)*t2(3)
        normal(3) = t1(1)*t2(2) - t1(2)*t2(1)
      else
        projpt    = 0
      endif

c     Format

3000  format(/' *ERROR* PRJ3DL: ',a)

      end
