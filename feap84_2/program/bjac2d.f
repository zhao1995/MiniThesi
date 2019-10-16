c$Id:$
      subroutine bjac2d ( rst , xl, ndm, shp, detj )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute jacobian determinant and shape functions
c               with natural derivatives for an 4-node quadrilateral.

c      Inputs:
c         rst(2)    - Natural coordinate location
c         xl(ndm,*) - Array of element coordinates
c         ndm       - Space dimension of mesh

c      Outputs:
c         shp(2,4)  - Shape function & derivative w/r natural coords.
c         detj      - Determinant of jacobian determinant
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm, i, j
      real*8    detj, xii, eti
      real*8    rst(2), xl(ndm,*), shp(3,4), xs(2,2), xi(4), eta(4)

      save

      data      xi /-0.5d0, 0.5d0, 0.5d0,-0.5d0/
      data      eta/-0.5d0,-0.5d0, 0.5d0, 0.5d0/

c     Compute shape functions and derivatives

      do i = 1,4
        xii      = 0.5d0 +  xi(i)*rst(1)
        eti      = 0.5d0 + eta(i)*rst(2)
        shp(1,i) = eti*xi(i)
        shp(2,i) = xii*eta(i)
        shp(3,i) = xii*eti
      end do ! i

c     Compute jacobian matrix

      do j = 1,2
        do i = 1,2
          xs(i,j) = xl(i,1)*shp(j,1) + xl(i,2)*shp(j,2)
     &            + xl(i,3)*shp(j,3) + xl(i,4)*shp(j,4)
        end do ! i
      end do ! j

c     Compute jacobian determinant

      detj = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)

      end
