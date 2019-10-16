c$Id:$
      subroutine shp3ds(ss, shp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 20-node serendipity shape functions

c      Inputs:
c         ss(*)    - Gauss points

c      Outputs:
c         shp(4,*) - Quadratic shape functions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    i, ix(4), iy(4), iz(4)
      real*8     ss(3), shp(4,*)
      real*8     xin(8),etn(8),ztn(8), m1(4),m2(4)
      real*8     xi0,xi1, et0,et1, zt0,zt1

      save

c     Ordering by vertex-edge

      data       ix /  9, 11, 13, 15 /
      data       iy / 12, 10, 16, 14 /
      data       iz / 17, 18, 20, 19 /

      data       xin / -1.d0,  1.d0,  1.d0, -1.d0,
     &                 -1.d0,  1.d0,  1.d0, -1.d0 /
      data       etn / -1.d0, -1.d0,  1.d0,  1.d0,
     &                 -1.d0, -1.d0,  1.d0,  1.d0 /
      data       ztn / -1.d0, -1.d0, -1.d0, -1.d0,
     &                  1.d0,  1.d0,  1.d0,  1.d0 /

      data       m1  / -1.d0,  1.d0, -1.d0,  1.d0 /
      data       m2  / -1.d0, -1.d0,  1.d0,  1.d0 /

c     Set 8-node corner shape functions

      do i = 1,8
        xi0 = xin(i)*ss(1)
        et0 = etn(i)*ss(2)
        zt0 = ztn(i)*ss(3)
        xi1 = 0.5d0 + 0.5d0*xi0
        et1 = 0.5d0 + 0.5d0*et0
        zt1 = 0.5d0 + 0.5d0*zt0
        shp(1,i) = et1*zt1*(ss(1) + 0.5d0*xin(i)*(et0 + zt0 - 1.0d0))
        shp(2,i) = xi1*zt1*(ss(2) + 0.5d0*etn(i)*(zt0 + xi0 - 1.0d0))
        shp(3,i) = xi1*et1*(ss(3) + 0.5d0*ztn(i)*(xi0 + et0 - 1.0d0))
        shp(4,i) = xi1*et1*zt1*(xi0 + et0 + zt0 - 2.d0)
      end do ! i

c     Mid-edge shape functions

      do i = 1,4
        xi1          = (1.d0 - ss(1)**2)*0.25d0
        et1          =  1.d0 + m1(i)*ss(2)
        zt1          =  1.d0 + m2(i)*ss(3)
        shp(1,ix(i)) = -0.5d0*ss(1)*et1*zt1
        shp(2,ix(i)) =  xi1*m1(i)*zt1
        shp(3,ix(i)) =  xi1*et1*m2(i)
        shp(4,ix(i)) =  xi1*et1*zt1

        xi1          =  1.d0 + m1(i)*ss(1)
        et1          = (1.d0 - ss(2)**2)*0.25d0
        shp(1,iy(i)) =  m1(i)*et1*zt1
        shp(2,iy(i)) = -0.5d0*xi1*ss(2)*zt1
        shp(3,iy(i)) =  xi1*et1*m2(i)
        shp(4,iy(i)) =  xi1*et1*zt1


        et1          =  1.d0 + m2(i)*ss(2)
        zt1          = (1.d0 - ss(3)**2)*0.25d0
        shp(1,iz(i)) =  m1(i)*et1*zt1
        shp(2,iz(i)) =  xi1*m2(i)*zt1
        shp(3,iz(i)) = -0.5d0*xi1*et1*ss(3)
        shp(4,iz(i)) =  xi1*et1*zt1
      end do ! i

      end
