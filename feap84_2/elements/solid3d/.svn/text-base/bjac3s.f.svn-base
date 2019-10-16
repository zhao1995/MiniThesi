c$Id:$
      subroutine bjac3s ( rst , xl, ndm, shp, detj )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute jacobian determinant of an 8-node brick

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm, i, j, k
      real*8    rst(3), xl(ndm,*), shp(4,8), xs(3,3)
      real*8    xi(8), eta(8), zta(8), detj, xii, eti, zti

      data xi /-0.5d0, 0.5d0, 0.5d0,-0.5d0,-0.5d0, 0.5d0, 0.5d0,-0.5d0/
      data eta/-0.5d0,-0.5d0, 0.5d0, 0.5d0,-0.5d0,-0.5d0, 0.5d0, 0.5d0/
      data zta/-0.5d0,-0.5d0,-0.5d0,-0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0/

c     Compute shape functions and their derivatives

      do i = 1,8
        xii = 0.5 +  xi(i)*rst(1)
        eti = 0.5 + eta(i)*rst(2)
        zti = 0.5 + zta(i)*rst(3)
        shp(1,i) =  xi(i)*eti*zti
        shp(2,i) = eta(i)*xii*zti
        shp(3,i) = zta(i)*xii*eti
        shp(4,i) =  xii*eti*zti
      end do ! i

c     Compute jacobian matrix

      do j = 1,3
        do i = 1,3
          xs(i,j) = 0.0d0
          do k = 1,8
            xs(i,j) = xs(i,j) + xl(i,k)*shp(j,k)
          end do ! k
        end do ! i
      end do ! j

c     Compute jacobian determinant

      detj = xs(1,1)*(xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2))
     &     + xs(1,2)*(xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3))
     &     + xs(1,3)*(xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1))

      end
