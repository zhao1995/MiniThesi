c$Id:$
      subroutine jacb3m(ss,xsj,xl,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Jacobian determinant of isoparametric map

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm, j
      real*8    xsj, ap1,am1,ap2,am2,ap3,am3, ad1,ad2,ad3
      real*8    ss(3),z(4,3),xl(ndm,8),xs(3,3)

c     Compute shape functions and their natural coord. derivatives

      ap1 = 1.0d0 + ss(1)
      am1 = 1.0d0 - ss(1)
      ap2 = 1.0d0 + ss(2)
      am2 = 1.0d0 - ss(2)
      ap3 = 1.0d0 + ss(3)
      am3 = 1.0d0 - ss(3)

c     Compute for ( - , - ) values

      z(1,3)  = am1*am2
      z(1,1)  = am2*am3
      z(1,2)  = am1*am3

c     Compute for ( + , + ) values

      z(2,3)  = ap1*ap2
      z(2,1)  = ap2*ap3
      z(2,2)  = ap1*ap3

c     Compute for ( - , + ) values

      z(3,3)  = am1*ap2
      z(3,1)  = am2*ap3
      z(3,2)  = am1*ap3

c     Compute for ( + , - ) values

      z(4,3)  = ap1*am2
      z(4,1)  = ap2*am3
      z(4,2)  = ap1*am3

c     Compute jacobian transformation

      do j = 1,3
        xs(1,j) =((xl(j,2) - xl(j,1))*z(1,1)
     &          + (xl(j,7) - xl(j,8))*z(2,1)
     &          + (xl(j,6) - xl(j,5))*z(3,1)
     &          + (xl(j,3) - xl(j,4))*z(4,1))*0.125d0
        xs(2,j) =((xl(j,4) - xl(j,1))*z(1,2)
     &          + (xl(j,7) - xl(j,6))*z(2,2)
     &          + (xl(j,8) - xl(j,5))*z(3,2)
     &          + (xl(j,3) - xl(j,2))*z(4,2))*0.125d0
        xs(3,j) =((xl(j,5) - xl(j,1))*z(1,3)
     &          + (xl(j,7) - xl(j,3))*z(2,3)
     &          + (xl(j,8) - xl(j,4))*z(3,3)
     &          + (xl(j,6) - xl(j,2))*z(4,3))*0.125d0
      end do ! j

c     Compute adjoint to jacobian

      ad1 = xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2)
      ad2 = xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3)
      ad3 = xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1)

c     Compute determinant of jacobian

      xsj  = xs(1,1)*ad1+xs(1,2)*ad2+xs(1,3)*ad3

      end
