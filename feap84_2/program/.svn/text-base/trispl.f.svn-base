c$Id:$
      subroutine trispl(el,xl,ndm, xsj,shp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Triangular 6-node B-spline shape function routine

c      Inputs:
c         el(3)     - Area coordinates for point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh

c      Outputs:
c         xsj       - Jacobian determinant at point
c         shp(3,*)  - Shape functions and derivatives at point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   ndm, i,j
      real*8    xsj,xsjr, x1,x2,x3, y1,y2,y3
      real*8    el1x,el2x,el3x, el1y,el2y,el3y
      real*8    el(3),xl(ndm,*), shp(3,*), xb(3,6)

      save

c     Set up adjusted co-ordinates for B-Splines

      do j = 1,ndm
        do i = 1,3
          xb(j,i) = xl(j,i)
        end do ! i
        xb(j,4) = 2.0d0*xl(j,4) - 0.5d0*(xl(j,1) + xl(j,2))
        xb(j,5) = 2.0d0*xl(j,5) - 0.5d0*(xl(j,2) + xl(j,3))
        xb(j,6) = 2.0d0*xl(j,6) - 0.5d0*(xl(j,3) + xl(j,1))
      end do ! j

c     Form Jacobian terms

      x1 = (xb(1,1)*el(1) + xb(1,4)*el(2) + xb(1,6)*el(3))*2.0d0
      x2 = (xb(1,2)*el(2) + xb(1,5)*el(3) + xb(1,4)*el(1))*2.0d0
      x3 = (xb(1,3)*el(3) + xb(1,6)*el(1) + xb(1,5)*el(2))*2.0d0

      y1 = (xb(2,1)*el(1) + xb(2,4)*el(2) + xb(2,6)*el(3))*2.0d0
      y2 = (xb(2,2)*el(2) + xb(2,5)*el(3) + xb(2,4)*el(1))*2.0d0
      y3 = (xb(2,3)*el(3) + xb(2,6)*el(1) + xb(2,5)*el(2))*2.0d0

      xsj  = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
      xsjr = 1.0d0
      if(xsj.ne.0.0d0) xsjr = 1.0d0/xsj
      xsj  = 0.5d0*xsj

c     Specify shape functions and their derivatives

      el1x = (y2-y3)*xsjr
      el1y = (x3-x2)*xsjr

      el2x = (y3-y1)*xsjr
      el2y = (x1-x3)*xsjr

      el3x = (y1-y2)*xsjr
      el3y = (x2-x1)*xsjr

c     Global shape functions and derivatives

      shp(1,1) = 2.0d0*el(1)*el1x
      shp(2,1) = 2.0d0*el(1)*el1y
      shp(3,1) = el(1)**2

      shp(1,2) = 2.0d0*el(2)*el2x
      shp(2,2) = 2.0d0*el(2)*el2y
      shp(3,2) = el(2)**2

      shp(1,3) = 2.0d0*el(3)*el3x
      shp(2,3) = 2.0d0*el(3)*el3y
      shp(3,3) = el(3)**2

      shp(1,4) = 2.0d0*(el(2)*el1x + el(1)*el2x)
      shp(2,4) = 2.0d0*(el(2)*el1y + el(1)*el2y)
      shp(3,4) = el(1)*el(2)*2.0d0

      shp(1,5) = 2.0d0*(el(3)*el2x + el(2)*el3x)
      shp(2,5) = 2.0d0*(el(3)*el2y + el(2)*el3y)
      shp(3,5) = el(2)*el(3)*2.0d0

      shp(1,6) = 2.0d0*(el(1)*el3x + el(3)*el1x)
      shp(2,6) = 2.0d0*(el(1)*el3y + el(3)*el1y)
      shp(3,6) = el(3)*el(1)*2.0d0

      end
