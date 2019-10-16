c$Id:$
      subroutine hsizend(xl, ndm,nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    17/03/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute maximum and minimum element size

c      Inputs:
c         xl(ndm,nel) - Element nodal coordinates
c         ndm         - Spatial dimension of mesh
c         nel         - Number of nodes on element

c      Outputs:
c         hsize(2)    - Element min/max size
c                       1 = minimum; 2 = maximum
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'qudshp.h'

      integer    ndm,nel, m,n,i
      real*8     xl(ndm,nel), dh

      do n = 1,nel-1
        do m = n+1,nel
          dh = (xl(1,n) - xl(1,m))**2
          do i = 2,ndm
            dh = dh + (xl(i,n) - xl(i,m))**2
          end do ! i
          dh = sqrt(dh)
          if(hsize(1).eq.0.0d0) then
            hsize(1) = dh
          else
            hsize(1) = min(hsize(1),dh)
          endif
          hsize(2) = max(hsize(2),dh)
        end do ! m
      end do ! n

      end
