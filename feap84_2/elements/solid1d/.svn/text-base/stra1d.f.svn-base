c$Id:$
      subroutine stra1d(d,xl,ul,tl,shp,ndf,ndm,nel,xx,ta,eps)

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
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),tl(*),shp(2,*)
      real*8    eps(9,3), xx,ta

      save

c     Compute strains and coordinates

      do j = 1,6
        eps(j,1) = 0.0d0
        eps(j,3) = 0.0d0
      end do ! j
      xx = 0.0d0
      ta = -d(9)
      do j = 1,nel
        xx = xx + shp(2,j)*xl(1,j)
        ta = ta + shp(2,j)*tl(j)
        eps(1,1) = eps(1,1) + shp(1,j)*ul(1,j,1)
        eps(3,1) = eps(3,1) + shp(2,j)*ul(1,j,1)
        eps(1,3) = eps(1,1) + shp(1,j)*ul(1,j,2)
        eps(3,3) = eps(3,1) + shp(2,j)*ul(1,j,2)
      end do ! j

c     Set reference and current coords

      xref(1) = xx
      xref(2) = 0.0d0
      xref(3) = 0.0d0
      xcur(1) = xx + eps(3,1)
      xcur(2) = 0.0d0
      xcur(3) = 0.0d0

c     Compute enhanced strains

      if(etype.eq.3) then
        do j = 1,2
          eps(1,1) = eps(1,1) + shpi(1,j)*ui(j,1)
          eps(1,3) = eps(1,3) + shpi(1,j)*ui(j,2)
        end do ! j
      endif

c     Strain at t_n

      eps(1,2) = eps(1,1) - eps(1,3)
      eps(3,2) = eps(3,1) - eps(3,3)

c     Set 3-strain (thickness/hoop)

      if(stype.eq.3) then
        eps(3,1) = eps(3,1)/xx
        eps(3,2) = eps(3,2)/xx
        eps(3,3) = eps(3,3)/xx
      else
        eps(3,1) = 0.0d0
        eps(3,2) = 0.0d0
        eps(3,3) = 0.0d0
      endif

      end
