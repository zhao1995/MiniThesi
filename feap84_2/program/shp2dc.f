c$Id:$
      subroutine shp2dc(ss,xl, shp, xsj, ord, flg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Compress constants for xi*s2 ans xi*s9           01/01/2014
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Shape function routine for cubic (16-node) elements
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

      include   'pconstant.h'
      include   'sdata.h'

      logical    flg
      integer    ord, i,j,k, xi1(16),xi2(16)
      real*8     ss(2),xl(ndm,*),xds(2,2), shp(3,16), xsj
      real*8     n1(4),dn1(4), n2(4),dn2(4)
      real*8     mn1, xi1s2,xi2s2,xi1s9,xi2s9, sq1,sq2

      save

      data       xi1/1,2,2,1,3,4,2,2,4,3,1,1,3,4,4,3/
      data       xi2/1,1,2,2,1,1,3,4,2,2,4,3,3,3,4,4/

c     Do Shape functions

      sq1    = ss(1)*ss(1)
      sq2    = ss(2)*ss(2)
      xi1s9  = 9.d0*sq1 - 1.d0
      xi2s9  = 9.d0*sq2 - 1.d0
      xi1s2  = 1.d0 - sq1
      xi2s2  = 1.d0 - sq2

      n1(1)  = (1.d0 - ss(1))*xi1s9*0.0625d0
      n1(2)  = (1.d0 + ss(1))*xi1s9*0.0625d0
      n1(3)  = (one3 - ss(1))*xi1s2*1.6875d0
      n1(4)  = (one3 + ss(1))*xi1s2*1.6875d0

      n2(1)  = (1.d0 - ss(2))*xi2s9*0.0625d0
      n2(2)  = (1.d0 + ss(2))*xi2s9*0.0625d0
      n2(3)  = (one3 - ss(2))*xi2s2*1.6875d0
      n2(4)  = (one3 + ss(2))*xi2s2*1.6875d0

      dn1(1) = (  1.d0 + (18.d0 - 27.d0*ss(1))*ss(1))*0.0625d0
      dn1(2) = ( -1.d0 + (18.d0 + 27.d0*ss(1))*ss(1))*0.0625d0
      dn1(3) = (-27.d0 - (18.d0 - 81.d0*ss(1))*ss(1))*0.0625d0
      dn1(4) = ( 27.d0 - (18.d0 + 81.d0*ss(1))*ss(1))*0.0625d0

      dn2(1) = (  1.d0 + (18.d0 - 27.d0*ss(2))*ss(2))*0.0625d0
      dn2(2) = ( -1.d0 + (18.d0 + 27.d0*ss(2))*ss(2))*0.0625d0
      dn2(3) = (-27.d0 - (18.d0 - 81.d0*ss(2))*ss(2))*0.0625d0
      dn2(4) = ( 27.d0 - (18.d0 + 81.d0*ss(2))*ss(2))*0.0625d0

      do k = 1,16
        shp(3,k) = n1(xi1(k))*n2(xi2(k))
      end do ! k

c     Do first derivatives

      if(ord.ge.1) then

c       Local derivatives

        do k = 1,16
          shp(1,k) = dn1(xi1(k))* n2(xi2(k))
          shp(2,k) =  n1(xi1(k))*dn2(xi2(k))
        end do ! k

c       Jacobian matrix

        do j = 1,2
          do i = 1,2
            xds(i,j) = 0.0d0
            do k = 1,16
              xds(i,j) = xds(i,j) + xl(i,k)*shp(j,k)
            end do ! k
          end do ! i
        end do ! j
        xsj = xds(1,1)*xds(2,2) - xds(1,2)*xds(2,1)

c       Global derivatives

        if(.not.flg) then
          do k = 1,16
            mn1      = ( xds(2,2)*shp(1,k) - xds(2,1)*shp(2,k))/xsj
            shp(2,k) = (-xds(1,2)*shp(1,k) + xds(1,1)*shp(2,k))/xsj
            shp(1,k) = mn1
          end do ! k
        endif

      endif

      end
