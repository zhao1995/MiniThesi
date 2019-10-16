c$Id:$
      subroutine fbody3d(d,xl, r, ndm,ndf, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add fixed end moments                            09/03/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute body force loads
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'elbody.h'
      include   'eldata.h'
      include   'pconstant.h'

      integer    ndm,ndf,isw, ii
      real*8     d(*),xl(ndm,*), r(ndf,*), body(3), le

c     Set body loading factors

      if(isw.eq.15) then
        do ii = 1,3
          body(ii) = bodyf(ii)
        end do ! ii
      else
        call sbodyf(d, body)
      endif

c     Add body forces

      le       = sqrt((xl(1,2)-xl(1,1))**2 + (xl(2,2)-xl(2,1))**2
     &              + (xl(3,2)-xl(3,1))**2)*0.5d0
      do ii = 1,3
        r(ii,1)   = r(ii,1) + body(ii)*le
        r(ii,2)   = r(ii,2) + body(ii)*le
      end do ! ii

c     Check for cubic displacements and compute fixed end moments

      if(d(79).gt.0.0d0) then

c       x-body loading

        if(body(1).ne.0.0d0) then
          le = xl(2,2) - xl(2,1)
          le = le*abs(le)*one12
          r(6,1) = r(6,1) - body(1)*le
          r(6,2) = r(6,2) + body(1)*le
          le = xl(3,2) - xl(3,1)
          le = le*abs(le)*one12
          r(5,1) = r(5,1) + body(1)*le
          r(5,2) = r(5,2) - body(1)*le
        endif

c       y-body loading

        if(body(2).ne.0.0d0) then
          le = xl(3,2) - xl(3,1)
          le = le*abs(le)*one12
          r(4,1) = r(4,1) - body(2)*le
          r(4,2) = r(4,2) + body(2)*le
          le = xl(1,2) - xl(1,1)
          le = le*abs(le)*one12
          r(6,1) = r(6,1) + body(2)*le
          r(6,2) = r(6,2) - body(2)*le
        endif

c       z-body loading

        if(body(3).ne.0.0d0) then
          le = xl(1,2) - xl(1,1)
          le = le*abs(le)*one12
          r(5,1) = r(5,1) - body(3)*le
          r(5,2) = r(5,2) + body(3)*le
          le = xl(2,2) - xl(2,1)
          le = le*abs(le)*one12
          r(4,1) = r(4,1) + body(3)*le
          r(4,2) = r(4,2) - body(3)*le
        endif

      endif ! End fixed end moment

      end
