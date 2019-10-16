c$Id:$
      subroutine xcor3d(ss,xl,ixl, x)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute shape functions and coordinates for each point

c      Inputs:
c         ss(3)   - Natural coordinates for point
c         xl(3,*) - Nodal coordinates for brick
c         ixl(*)  - List of active nodes on brick

c      Outputs:
c         x(3)    - Cartesian coordinates for point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   j,l, ixl(27)
      real*8    ss(3),shp(27),xl(3,27),x(3)

      save

      call gshp3(ss, ixl, shp)

      do j = 1,3
        x(j) = 0.0d0
      end do ! j

      do l = 1,27
        if(ixl(l).ne.0) then
          do j = 1,3
            x(j) = x(j) + shp(l)*xl(j,l)
          end do ! j
        endif
      end do ! l

      end
