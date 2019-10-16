c$Id:$
      subroutine shp2ds(ss,xl, shp, xsj, ord, flg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Shape function routine for cubic (12-node) elements
c     Inputs:
c       ss(*)    - Gauss point
c       xl(*)    - Element coordinates
c       ord      - Order to generate:
c       flg      - .false. for global derivatives

c     Outputs:
c       shp(*)   - Shape functions    (ord.ge.0) and
c                  first derivatives  (ord.ge.1)
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'sdata.h'

      logical    flg
      integer    ord, i,j,k
      real*8     xds(2,2)
      real*8     ss(2),xl(ndm,*),shp(3,16), xsj
      real*8     n1(2),dn1(2), n2(2),dn2(2), xi1a(4), xi2a(4)
      real*8     mn1, xi1s2,xi2s2

      data       xi1a / -1.d0, 1.d0,1.d0,-1.d0/
      data       xi2a / -1.d0,-1.d0,1.d0, 1.d0/

c     Do Shape functions

      xi1s2  = 1.d0 - ss(1)*ss(1)
      xi2s2  = 1.d0 - ss(2)*ss(2)

      n1(1)  = 9.d0*xi1s2*(1.d0 - 3.d0*ss(1))*0.0625d0
      n1(2)  = 9.d0*xi1s2*(1.d0 + 3.d0*ss(1))*0.0625d0

      n2(1)  = 9.d0*xi2s2*(1.d0 - 3.d0*ss(2))*0.0625d0
      n2(2)  = 9.d0*xi2s2*(1.d0 + 3.d0*ss(2))*0.0625d0

      dn1(1) = (-27.d0 - (18.d0 - 81.d0*ss(1))*ss(1))*0.0625d0
      dn1(2) = ( 27.d0 - (18.d0 + 81.d0*ss(1))*ss(1))*0.0625d0

      dn2(1) = (-27.d0 - (18.d0 - 81.d0*ss(2))*ss(2))*0.0625d0
      dn2(2) = ( 27.d0 - (18.d0 + 81.d0*ss(2))*ss(2))*0.0625d0

      do i = 1,4
        shp(1,i) = 0.25d0*xi1a(i)*(1.d0 + xi2a(i)*ss(2))
        shp(2,i) = 0.25d0*xi2a(i)*(1.d0 + xi1a(i)*ss(1))
        shp(3,i) = 0.25d0*(1.d0 + xi1a(i)*ss(1))*(1.d0 + xi2a(i)*ss(2))
      end do ! i

      xi1s2     = 0.5d0*(1.d0 - ss(1))
      xi2s2     = 0.5d0*(1.d0 - ss(2))

      shp(1,5)  = dn1(1)*xi2s2
      shp(2,5)  = -n1(1)*0.5d0
      shp(3,5)  =  n1(1)*xi2s2

      shp(1,6)  = dn1(2)*xi2s2
      shp(2,6)  = -n1(2)*0.5d0
      shp(3,6)  =  n1(2)*xi2s2

      shp(1,12) = -n2(1)*0.5d0
      shp(2,12) = dn2(1)*xi1s2
      shp(3,12) =  n2(1)*xi1s2

      shp(1,11) = -n2(2)*0.5d0
      shp(2,11) = dn2(2)*xi1s2
      shp(3,11) =  n2(2)*xi1s2

      xi1s2     = 0.5d0*(1.d0 + ss(1))
      xi2s2     = 0.5d0*(1.d0 + ss(2))

      shp(1,10)  = dn1(1)*xi2s2
      shp(2,10)  =  n1(1)*0.5d0
      shp(3,10)  =  n1(1)*xi2s2

      shp(1, 9)  = dn1(2)*xi2s2
      shp(2, 9)  =  n1(2)*0.5d0
      shp(3, 9)  =  n1(2)*xi2s2

      shp(1, 7)  =  n2(1)*0.5d0
      shp(2, 7)  = dn2(1)*xi1s2
      shp(3, 7)  =  n2(1)*xi1s2

      shp(1, 8)  =  n2(2)*0.5d0
      shp(2, 8)  = dn2(2)*xi1s2
      shp(3, 8)  =  n2(2)*xi1s2

      do i = 1,3
        shp(i, 1)  = shp(i,1) - (2.d0*(shp(i, 5) + shp(i,12))
     &                               + shp(i, 6) + shp(i,11))
        shp(i, 2)  = shp(i,2) - (2.d0*(shp(i, 6) + shp(i, 7))
     &                               + shp(i, 5) + shp(i, 8))
        shp(i, 3)  = shp(i,3) - (2.d0*(shp(i, 8) + shp(i, 9))
     &                               + shp(i, 7) + shp(i,10))
        shp(i, 4)  = shp(i,4) - (2.d0*(shp(i,10) + shp(i,11))
     &                               + shp(i, 9) + shp(i,12))
      end do ! i

      if(ord.ge.1) then

c       Jacobian matrix

        do j = 1,2
          do i = 1,2
            xds(i,j) = 0.0d0
            do k = 1,12
              xds(i,j) = xds(i,j) + xl(i,k)*shp(j,k)
            end do ! k
          end do ! i
        end do ! j
        xsj = xds(1,1)*xds(2,2) - xds(1,2)*xds(2,1)

c       Global derivatives

        if(.not.flg) then
          do k = 1,12
            mn1      = ( xds(2,2)*shp(1,k) - xds(2,1)*shp(2,k))/xsj
            shp(2,k) = (-xds(1,2)*shp(1,k) + xds(1,1)*shp(2,k))/xsj
            shp(1,k) = mn1
          end do ! k
        endif

      endif

      end
