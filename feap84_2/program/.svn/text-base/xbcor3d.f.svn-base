c$Id:$
      subroutine xbcor3d(ss,xl, x)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise shape function numbering to standare      10/11/2008
c          27 node Lagrange element
c       2. Restore option for old numbering on flag oldfl   28/01/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute shape functions and coordinates for each point

c      Inputs:
c         ss(3)   - Natural coordinates for point
c         xl(3,*) - Nodal coordinates for brick

c      Outputs:
c         x(3)    - Cartesian coordinates for point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'corset.h'

      integer   j,l, ix(27),iy(27),iz(27), ox(27),oy(27),oz(27)
      real*8    ss(3),xl(3,27),x(3), lshp(3,3),shp

      save

      data      ix/1,3,3,1, 1,3,3,1, 2,3,2,1, 2,3,2,1, 1,3,3,1,
     &             1,3, 2,2, 2,2, 2/
      data      iy/1,1,3,3, 1,1,3,3, 1,2,3,2, 1,2,3,2, 1,1,3,3,
     &             2,2, 3,1, 2,2, 2/
      data      iz/1,1,1,1, 3,3,3,3, 1,1,1,1, 3,3,3,3, 2,2,2,2,
     &             2,2, 2,2, 1,3, 2/

      data      ox/1,3,3,1, 1,3,3,1, 1,3,3,1, 2,3,2,1,2, 2,3,2,1,2,
     &             2,3,2,1,2/
      data      oy/1,1,3,3, 1,1,3,3, 1,1,3,3, 1,2,3,2,2, 1,2,3,2,2,
     &             1,2,3,2,2/
      data      oz/1,1,1,1, 3,3,3,3, 2,2,2,2, 1,1,1,1,1, 3,3,3,3,3,
     &             2,2,2,2,2/


      do j = 1,3
        lshp(1,j) = 0.5d0*ss(j)*(ss(j) - 1.d0)
        lshp(2,j) = (1.d0 - ss(j)*ss(j))
        lshp(3,j) = 0.5d0*ss(j)*(ss(j) + 1.d0)
      end do ! j

      do j = 1,3
        x(j) = 0.0d0
      end do ! j

      do l = 1,27
        if(oldfl) then
          shp = lshp(ox(l),1)*lshp(oy(l),2)*lshp(oz(l),3)
        else
          shp = lshp(ix(l),1)*lshp(iy(l),2)*lshp(iz(l),3)
        endif
        do j = 1,3
          x(j) = x(j) + shp*xl(j,l)
        end do ! j
      end do ! l

      end
