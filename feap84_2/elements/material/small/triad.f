c$Id:$
      subroutine triad(r)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    20/12/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute ON triad from two vectors

c      Inputs:
c         r(3,1)   - vector 1
c         r(3,2)   - vector 2

c      Outputs:
c         r(3,1)   - unit length
c         r(3,2)   - projection of b - b.a a, unit length
c         r(3,3)   - a x b
c         r(3,3) - rotation matrix from a,b,c to e1,e2,e3

c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      real*8     r(3,3)
      real*8     nv

c     Normalize r_1
      nv = 1.0d0/sqrt(r(1,1)*r(1,1)+r(2,1)*r(2,1)+r(3,1)*r(3,1))
      r(1,1) = r(1,1)*nv
      r(2,1) = r(2,1)*nv
      r(3,1) = r(3,1)*nv

c     Project r_2 to r_1-perp
      nv = r(1,1)*r(1,2)+r(2,1)*r(2,2)+r(3,1)*r(3,2)
      r(1,2) = r(1,2) - nv*r(1,1)
      r(2,2) = r(2,2) - nv*r(2,1)
      r(3,2) = r(3,2) - nv*r(3,1)

c     Normalize r_2
      nv = 1.0d0/sqrt(r(1,2)*r(1,2)+r(2,2)*r(2,2)+r(3,2)*r(3,2))
      r(1,2) = r(1,2)*nv
      r(2,2) = r(2,2)*nv
      r(3,2) = r(3,2)*nv

c     get r_3 = r_1 x r_2
      r(1,3) = r(2,1)*r(3,2)-r(3,1)*r(2,2)
      r(2,3) = r(3,1)*r(1,2)-r(1,1)*r(3,2)
      r(3,3) = r(1,1)*r(2,2)-r(2,1)*r(1,2)

c     Normalize r_3
      nv = 1.0d0/sqrt(r(1,3)*r(1,3)+r(2,3)*r(2,3)+r(3,3)*r(3,3))
      r(1,3) = r(1,3)*nv
      r(2,3) = r(2,3)*nv
      r(3,3) = r(3,3)*nv

      end
