c$Id:$
      subroutine tetshph( xi, xl, ndm, order, xsj, shp, hierfl )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Program shape functions for L_1=1-L_2-L_3-L_4    13/08/2007
c          Divide jacobian by 6.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute 3-d tetrahedral element shape
c               functions and their derivatives w/r x,y,z

c      Inputs:
c         xi(4)     - Natural volume coordinates of point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c         order     - Interpolation order: 1=linear; 2=quadratic
c         hierfl    - Use hierarchic interpolation if true

c      Outputs:
c         xsj       - Jacobian determinant at point
c         shp(4,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dx
c                     shp(2,i) = dN_i/dy
c                     shp(3,i) = dN_i/dz
c                     shp(4,i) =  N_i
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'eldata.h'
      include 'iofile.h'
      include 'pconstant.h'

      logical  hierfl
      integer  ndm, order, i, j, k
      real*8   xsj, detr, xi(4), xl(ndm,*), shp(4,*), xs(3,3),xsi(3,3)
      real*8   te(3)

      save

c     Linear shape functions and their derivatives

      if(order.eq.1) then
        shp(1, 1) =-1.d0
        shp(1, 2) = 1.d0
        shp(1, 3) = 0.d0
        shp(1, 4) = 1.d0

        shp(2, 1) =-1.d0
        shp(2, 2) = 0.d0
        shp(2, 3) = 1.d0
        shp(2, 4) = 0.d0

        shp(3, 1) =-1.d0
        shp(3, 2) = 0.d0
        shp(3, 3) = 0.d0
        shp(3, 4) = 1.d0

        shp(4, 1) = xi(1)
        shp(4, 2) = xi(2)
        shp(4, 3) = xi(3)
        shp(4, 4) = xi(4)

c     Quadratic shape functions and derivatives

      elseif(order.eq.2) then

c       Shape functions for vertex nodes

        if(hierfl) then
          shp(1, 1) =-1.d0
          shp(1, 2) = 1.d0
          shp(1, 3) = 0.d0
          shp(1, 4) = 0.d0

          shp(2, 1) =-1.d0
          shp(2, 2) = 0.d0
          shp(2, 3) = 1.d0
          shp(2, 4) = 0.d0

          shp(3, 1) =-1.d0
          shp(3, 2) = 0.d0
          shp(3, 3) = 0.d0
          shp(3, 4) = 1.d0

          shp(4, 1) = xi(1)
          shp(4, 2) = xi(2)
          shp(4, 3) = xi(3)
          shp(4, 4) = xi(4)
        else
          shp(1, 1) =-4.d0*xi(1) + 1.d0
          shp(1, 2) = 4.d0*xi(2) - 1.d0
          shp(1, 3) = 0.d0
          shp(1, 4) = 0.d0

          shp(2, 1) =-4.d0*xi(1) + 1.d0
          shp(2, 2) = 0.d0
          shp(2, 3) = 4.d0*xi(3) - 1.d0
          shp(2, 4) = 0.d0

          shp(3, 1) =-4.d0*xi(1) + 1.d0
          shp(3, 2) = 0.d0
          shp(3, 3) = 0.d0
          shp(3, 4) = 4.d0*xi(4) - 1.d0

          shp(4, 1) = xi(1)*(2.d0*xi(1) - 1.d0)
          shp(4, 2) = xi(2)*(2.d0*xi(2) - 1.d0)
          shp(4, 3) = xi(3)*(2.d0*xi(3) - 1.d0)
          shp(4, 4) = xi(4)*(2.d0*xi(4) - 1.d0)
        endif

c       Shape functions for mid-edge nodes

        shp(1, 5) = 4.d0*(xi(1) - xi(2))
        shp(1, 6) = 4.d0*xi(3)
        shp(1, 7) =-4.d0*xi(3)
        shp(1, 8) =-4.d0*xi(4)
        shp(1, 9) = 4.d0*xi(4)
        shp(1,10) = 0.d0

        shp(2, 5) =-4.d0*xi(2)
        shp(2, 6) = 4.d0*xi(2)
        shp(2, 7) = 4.d0*(xi(1) - xi(3))
        shp(2, 8) =-4.d0*xi(4)
        shp(2, 9) =-0.d0
        shp(2,10) = 4.d0*xi(4)

        shp(3, 5) =-4.d0*xi(2)
        shp(3, 6) = 0.d0
        shp(3, 7) =-4.d0*xi(3)
        shp(3, 8) = 4.d0*(xi(1) - xi(4))
        shp(3, 9) = 4.d0*xi(2)
        shp(3,10) = 4.d0*xi(3)

        shp(4, 5) = 4.d0*xi(1)*xi(2)
        shp(4, 6) = 4.d0*xi(2)*xi(3)
        shp(4, 7) = 4.d0*xi(3)*xi(1)
        shp(4, 8) = 4.d0*xi(1)*xi(4)
        shp(4, 9) = 4.d0*xi(2)*xi(4)
        shp(4,10) = 4.d0*xi(3)*xi(4)

        if(nel.eq.11) then
          shp(1,11) = 256.d0*(xi(1) - xi(2))*xi(3)*xi(4)
          shp(2,11) = 256.d0*(xi(1) - xi(3))*xi(2)*xi(4)
          shp(3,11) = 256.d0*(xi(1) - xi(4))*xi(2)*xi(3)
          shp(4,11) = 256.d0*xi(1)*xi(2)*xi(3)*xi(4)
        endif

c     Error - Higher than quadratic not coded

      else

        write(iow,2000) order
        write(  *,2000) order
        call plstop()

      endif

c     Compute jacobian matrix

      do i = 1,3
        do j = 1,3
          xs(j,i) = 0.0d0
          do k = 1,nel
            xs(j,i) = xs(j,i) + xl(j,k)*shp(i,k)
          end do ! k
        end do ! i
      end do ! j

c     Compute inverse of jacobian matrix

      xsi(1,1) = xs(2,2)*xs(3,3) - xs(3,2)*xs(2,3)
      xsi(1,2) = xs(3,2)*xs(1,3) - xs(1,2)*xs(3,3)
      xsi(1,3) = xs(1,2)*xs(2,3) - xs(2,2)*xs(1,3)

      xsi(2,1) = xs(2,3)*xs(3,1) - xs(3,3)*xs(2,1)
      xsi(2,2) = xs(3,3)*xs(1,1) - xs(1,3)*xs(3,1)
      xsi(2,3) = xs(1,3)*xs(2,1) - xs(2,3)*xs(1,1)

      xsi(3,1) = xs(2,1)*xs(3,2) - xs(3,1)*xs(2,2)
      xsi(3,2) = xs(3,1)*xs(1,2) - xs(1,1)*xs(3,2)
      xsi(3,3) = xs(1,1)*xs(2,2) - xs(2,1)*xs(1,2)

c     Compute jacobian determinant

      xsj = xs(1,1)*xsi(1,1) + xs(1,2)*xsi(2,1) + xs(1,3)*xsi(3,1)

      if(xsj.ne.0.0d0) then
        detr = 1.d0/xsj
        xsj  = xsj*one6
      else
        write(iow,*) ' TETSHPH: Determinant =',xsj
        detr = 1.d0
      endif

c     Compute jacobian inverse

      do j = 1,3
        do i = 1,3
          xs(i,j) = xsi(i,j)*detr
        end do ! i
      end do ! j

c     Compute shape function derivatives

      do k = 1,nel
        do i = 1,3
          te(i) = shp(1,k)*xs(1,i) + shp(2,k)*xs(2,i) + shp(3,k)*xs(3,i)
        end do ! i
        do i = 1,3
          shp(i,k) = te(i)
        end do ! i
      end do ! k

c     Format

2000  format(/' *ERROR* TETSHP not coded for order =',i4)

      end
