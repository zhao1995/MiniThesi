c$Id:$
      subroutine tri10(el,xl,ndm,xsj,shps)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct computation of jacobian and loop on shps 02/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Shape function routine for 10-node triangle

c      Inputs:
c         el(3)     - Area coordinates for point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh

c      Outputs:
c         xsj       - Jacobian determinant at point
c         shp(3,*)  - Shape functions and derivatives at point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    ndm
      real*8     el(3),xl(ndm,*), xsj, shps(3,*), dshp(3,10)
      real*8     xs(3,3)

      integer    i,j,k,n

      n = 3
      do i = 1,3
        j        = mod(i,3) + 1
        k        = mod(j,3) + 1

c       Vertex shape functions (1,2,3)

        shps(3,i) = 4.5d0*el(i)*el(i)*(el(i) - 1.d0) + el(i)
        dshp(i,i) = (13.5d0*el(i) - 9.0d0)*el(i) + 1.0d0
        dshp(j,i) = 0.0d0
        dshp(k,i) = 0.0d0

c       Side shape functions (4,5, 6,7, 8,9)

        n         = n + 1
        shps(3,n) = 4.5d0*el(i)*el(j)*(3.d0*el(i) - 1.0d0)
        dshp(i,n) = 4.5d0*el(j)*(6.d0*el(i) - 1.0d0)
        dshp(j,n) = 4.5d0*el(1)*(3.d0*el(i) - 1.0d0)
        dshp(k,n) = 0.0d0
        n         = n + 1
        shps(3,n) = 4.5d0*el(i)*el(j)*(3.d0*el(j) - 1.0d0)
        dshp(i,n) = 4.5d0*el(j)*(3.d0*el(j) - 1.0d0)
        dshp(j,n) = 4.5d0*el(i)*(6.d0*el(j) - 1.0d0)
        dshp(k,n) = 0.0d0
      end do ! i

c     Baricenter shape function

      shps(3,10) = 27.d0*el(1)*el(2)*el(3)
      dshp(1,10) = 27.d0*el(1)*el(2)
      dshp(2,10) = 27.d0*el(2)*el(3)
      dshp(3,10) = 27.d0*el(3)*el(1)

c     Compute the jacobian

      do j = 1,3
        do i = 1,2
          xs(i,j) = 0.0d0
          do k = 1,10
            xs(i,j) = xs(i,j) + xl(i,k)*dshp(j,k)
          end do ! k
        end do ! i
        xs(3,j) = 1.0d0
      end do ! j
      xsj = xs(1,1)*(xs(2,2) - xs(2,3))
     &    + xs(1,2)*(xs(2,3) - xs(2,1))
     &    + xs(1,3)*(xs(2,1) - xs(2,2))
      xsj = xsj*0.5d0

      call invert(xs,3,3)
      do k = 1,10
        do i = 1,2
          shps(i,k) = dshp(1,k)*xs(1,i)
     &              + dshp(2,k)*xs(2,i)
     &              + dshp(3,k)*xs(3,i)
        end do ! i
      end do ! k

      end
