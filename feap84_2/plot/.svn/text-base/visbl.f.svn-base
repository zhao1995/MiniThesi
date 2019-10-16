c$Id:$
      logical function visbl(xl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Find a visible face for quadrilateral facets

c      Inputs:
c         xl(3,*)   - Nodal coordinates for facet

c      Outputs:
c         visbl     - Flag, true if visible from view point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ppers.h'
      include  'pdatay.h'

      integer   k, jsym
      real*8    xl(3,4),v1(3),v2(3),v3(3), dot

      save

      do k = 1,3
        v1(k) = xl(k,3) - xl(k,1)
        v2(k) = xl(k,4) - xl(k,2)
      end do ! k

      call vecp(v1,v2,v3)

      do k = 1,3
        v1(k) = e(k) - (xl(k,1)+xl(k,2)+xl(k,3)+xl(k,4))*0.25d0
      end do ! k

      jsym  = max(1,lsym)

      visbl = (vsym(jsym)*dot(v3,v1,3)).gt.0.0d0

      end
