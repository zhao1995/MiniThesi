c$Id:$
      subroutine strn2d(d,xl,ul,tl,shp,ndf,ndm,nel,xx,yy,ta,eps)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'elcoor.h'
      include  'incshp.h'
      include  'pmod2d.h'

      integer   ndf,ndm,nel, j
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),tl(*),shp(3,*)
      real*8    eps(9,3), xx,yy,uu,vv(3),ww,ta

      save

c     Compute strains and coordinates

      do j = 1,6
        eps(j,1) = 0.0d0
        eps(j,3) = 0.0d0
      end do ! j
      xx = 0.0d0
      yy = 0.0d0
      uu = 0.0d0
      ww = 0.0d0
      ta = -d(9)
      do j = 1,nel
        xx = xx + shp(3,j)*xl(1,j)
        yy = yy + shp(3,j)*xl(2,j)
        uu = uu + shp(3,j)*ul(1,j,1)
        ww = ww + shp(3,j)*ul(2,j,1)
        ta = ta + shp(3,j)*tl(j)
        eps(1,1) = eps(1,1) + shp(1,j)*ul(1,j,1)
        eps(2,1) = eps(2,1) + shp(2,j)*ul(2,j,1)
        eps(3,1) = eps(3,1) + shp(3,j)*ul(1,j,1)
        eps(4,1) = eps(4,1) + shp(2,j)*ul(1,j,1) + shp(1,j)*ul(2,j,1)
        eps(1,3) = eps(1,3) + shp(1,j)*ul(1,j,2)
        eps(2,3) = eps(2,3) + shp(2,j)*ul(2,j,2)
        eps(3,3) = eps(3,3) + shp(3,j)*ul(1,j,2)
        eps(4,3) = eps(4,3) + shp(2,j)*ul(1,j,2) + shp(1,j)*ul(2,j,2)
      end do ! j

c     Set reference and current coords

      xref(1) = xx
      xref(2) = yy
      xref(3) = 0.0d0
      xcur(1) = xx + uu
      xcur(2) = yy + ww
      xcur(3) = 0.0d0

c     Axisymmetry with torsion

      if(stype.eq.8) then
        vv(1) = 0.0d0
        vv(3) = 0.0d0
        do j = 1,nel
          eps(5,1) = eps(5,1) + shp(2,j)*ul(3,j,1)
          eps(6,1) = eps(6,1) + shp(1,j)*ul(3,j,1)
          vv(1)    = vv(1)    + shp(3,j)*ul(3,j,1)

          eps(5,3) = eps(5,3) + shp(2,j)*ul(3,j,2)
          eps(6,3) = eps(6,3) + shp(1,j)*ul(3,j,2)
          vv(3)    = vv(3)    + shp(3,j)*ul(3,j,2)
        end do ! j
        eps(6,1) = eps(6,1) - vv(1)/xx
        eps(6,3) = eps(6,3) - vv(3)/xx
      endif

c     Compute enhanced strains (none for torsion)

      if(etype.eq.3) then
        do j = 1,2
          eps(1,1) = eps(1,1) + shpi(1,j)*ui(2*j-1,1)
          eps(2,1) = eps(2,1) + shpi(2,j)*ui(2*j,1)
          eps(4,1) = eps(4,1) + shpi(1,j)*ui(2*j,1)
     &                        + shpi(2,j)*ui(2*j-1,1)
          eps(1,3) = eps(1,3) + shpi(1,j)*ui(2*j-1,2)
          eps(2,3) = eps(2,3) + shpi(2,j)*ui(2*j,2)
          eps(4,3) = eps(4,3) + shpi(1,j)*ui(2*j,2)
     &                        + shpi(2,j)*ui(2*j-1,2)
        end do ! j
        eps(3,1) = eps(3,1) + shpi(3,3)*ui(5,1)
        eps(3,3) = eps(3,3) + shpi(3,3)*ui(5,2)
      endif

c     Strain at t_n

      do j = 1,6
        eps(j,2) = eps(j,1) - eps(j,3)
      end do ! j

c     Set 3-strain (thickness/hoop)

      if(stype.eq.3 .or. stype.eq.8) then
        eps(3,1) = eps(3,1)/xx
        eps(3,2) = eps(3,2)/xx
        eps(3,3) = eps(3,3)/xx
      else
        eps(3,1) = 0.0d0
        eps(3,2) = 0.0d0
        eps(3,3) = 0.0d0
      endif

      end
