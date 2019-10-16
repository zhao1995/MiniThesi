c$Id:$
      subroutine xiso3d(xi1,xi2,xi3,xl,xx)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form coordinates in 8-node block

c      Inputs:
c        xi1          - xi   - coordinate value
c        xi2          - eta  - coordinate value
c        xi3          - zeta - coordinate value
c        xl(3,*)      - Nodal coordinates of block

c      Outputs:
c        xx(3)        - Coordinate of point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i,j
      real*8   sh1m,sh2m,sh3m,sh1p,sh2p,sh3p
      real*8   xi1,xi2,xi3, shp(8),xl(3,8), xx(3)

      save

c     Constant parameters for 8-node shape functions

      sh1m = 0.25d0 - 0.25d0*xi1
      sh1p = 0.25d0 + 0.25d0*xi1
      sh2m = 0.50d0 - 0.50d0*xi2
      sh2p = 0.50d0 + 0.50d0*xi2
      sh3m = 0.50d0 - 0.50d0*xi3
      sh3p = 0.50d0 + 0.50d0*xi3

c     Form shape functions

      shp(1) = sh1m*sh2m*sh3m
      shp(2) = sh1p*sh2m*sh3m
      shp(3) = sh1p*sh2p*sh3m
      shp(4) = sh1m*sh2p*sh3m
      shp(5) = sh1m*sh2m*sh3p
      shp(6) = sh1p*sh2m*sh3p
      shp(7) = sh1p*sh2p*sh3p
      shp(8) = sh1m*sh2p*sh3p

c     Subtract 2 x isoparametric interpolation coordinate

      do j = 1,3
        do i = 1,8
          xx(j) = xx(j) - shp(i)*xl(j,i)
        end do ! i
      end do ! j

      end
