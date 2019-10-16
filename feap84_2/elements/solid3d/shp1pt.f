c$Id:$
      subroutine shp1pt(xsj1,shp1,xji,xl,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  One-point quadrature shape function routine

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm, i , i0, i1, i2, j0, j1, j2
      real*8    shp1(4,8),xl(ndm,8),xj(3,3),xji(3,3), xsj1, xsji

      do i = 1,3

        xj(1,i) = 0.125d0*(- xl(i,1) + xl(i,2) + xl(i,3) - xl(i,4)
     &                     - xl(i,5) + xl(i,6) + xl(i,7) - xl(i,8))
        xj(2,i) = 0.125d0*(- xl(i,1) - xl(i,2) + xl(i,3) + xl(i,4)
     &                     - xl(i,5) - xl(i,6) + xl(i,7) + xl(i,8))
        xj(3,i) = 0.125d0*(- xl(i,1) - xl(i,2) - xl(i,3) - xl(i,4)
     &                     + xl(i,5) + xl(i,6) + xl(i,7) + xl(i,8))

      end do ! i

c     Compute element jacobian adjoint at center

      do i0 = 1,3
        i1 = mod(i0,3) + 1
        i2 = mod(i1,3) + 1
        do j0 = 1,3
          j1 = mod(j0,3) + 1
          j2 = mod(j1,3) + 1

          xji(i0,j0) = xj(j1,i1)*xj(j2,i2) - xj(j1,i2)*xj(j2,i1)

        end do ! j0
      end do ! i0

c     Compute element jacobian determinant and its inverse

      xsj1 = xj(1,1)*xji(1,1) + xj(1,2)*xji(2,1) + xj(1,3)*xji(3,1)

      xsji = 0.125d0/xsj1

c     Compute element shape function derivs (current configuration)

      do i = 1,3

        shp1(i,1) = ( - xji(i,1) - xji(i,2) - xji(i,3) )*xsji
        shp1(i,2) = ( + xji(i,1) - xji(i,2) - xji(i,3) )*xsji
        shp1(i,3) = ( + xji(i,1) + xji(i,2) - xji(i,3) )*xsji
        shp1(i,4) = ( - xji(i,1) + xji(i,2) - xji(i,3) )*xsji

c       Other four by asymmetry

        shp1(i,5) = - shp1(i,3)
        shp1(i,6) = - shp1(i,4)
        shp1(i,7) = - shp1(i,1)
        shp1(i,8) = - shp1(i,2)

      end do ! i

c     Compute shape functions

      do i = 1,8
        shp1(4,i) = 0.125d0
      end do ! i

      end
