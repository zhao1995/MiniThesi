c$Id:$
      subroutine strn3d(d,xl,ul,th,shp,ndm,ndf,nel, eps,ta)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Three dimensional strain calculations

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'elcoor.h'

      integer   ndm,ndf,nel, i,j
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),th(*),shp(4,*)
      real*8    eps(9,*), ta

c     Compute stress and strain at point

      do j = 1,9
        eps(j,1) = 0.0d0
        eps(j,2) = 0.0d0
        eps(j,3) = 0.0d0
      end do ! j

      do i = 1,3
        xref(i) = 0.0d0
        xcur(i) = 0.0d0
      end do ! i

c     Compute temperature and coordinates

      ta = -d(9)
      do j = 1,nel
        ta     = ta         + shp(4,j)*th(j)
        do i = 1,3
          xref(i)  = xref(i)  + shp(4,j)*xl(i,j)
          xcur(i)  = xcur(i)  + shp(4,j)*ul(i,j,1)
        end do ! i
        eps(1,1) = eps(1,1) + shp(1,j)*ul(1,j,1)
        eps(2,1) = eps(2,1) + shp(2,j)*ul(2,j,1)
        eps(3,1) = eps(3,1) + shp(3,j)*ul(3,j,1)
        eps(4,1) = eps(4,1) + shp(1,j)*ul(2,j,1)
     &                      + shp(2,j)*ul(1,j,1)
        eps(5,1) = eps(5,1) + shp(2,j)*ul(3,j,1)
     &                      + shp(3,j)*ul(2,j,1)
        eps(6,1) = eps(6,1) + shp(3,j)*ul(1,j,1)
     &                      + shp(1,j)*ul(3,j,1)
        eps(1,3) = eps(1,3) + shp(1,j)*ul(1,j,2)
        eps(2,3) = eps(2,3) + shp(2,j)*ul(2,j,2)
        eps(3,3) = eps(3,3) + shp(3,j)*ul(3,j,2)
        eps(4,3) = eps(4,3) + shp(1,j)*ul(2,j,2) + shp(2,j)*ul(1,j,2)
        eps(5,3) = eps(5,3) + shp(2,j)*ul(3,j,2) + shp(3,j)*ul(2,j,2)
        eps(6,3) = eps(6,3) + shp(3,j)*ul(1,j,2) + shp(1,j)*ul(3,j,2)
      end do ! j
      do j = 1,6
        eps(j,2) = eps(j,1) - eps(j,3)
      end do ! j

      do i = 1,3
        xcur(i) = xcur(i) + xref(i)
      end do ! i

      end
