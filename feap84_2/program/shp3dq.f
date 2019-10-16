c$Id:$
      subroutine shp3dq(ss, shp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:
c         ss(*)    - Gauss points

c      Outputs:
c         shp(4,*) - Quadratic shape functions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    i, ir(27), is(27), it(27)
      real*8     nr(3),dr(3), ns(3),ds(3), nt(3),dt(3), ss(3), shp(4,*)

      save

c     Ordering by vertex-edge-face-interior

      data       ir /1,2,2,1, 1,2,2,1, 3,2,3,1,
     &               3,2,3,1, 1,2,2,1, 3,3,1,2,3,3, 3/
      data       is /1,1,2,2, 1,1,2,2, 1,3,2,3,
     &               1,3,2,3, 1,1,2,2, 3,3,3,3,1,2, 3/
      data       it /1,1,1,1, 2,2,2,2, 1,1,1,1,
     &               2,2,2,2, 3,3,3,3, 1,2,3,3,3,3, 3/

c     Set 1-d shape functions for each local direction

      nr(1) =  0.5d0 * ss(1) * (ss(1) - 1.d0)
      nr(2) =  0.5d0 * ss(1) * (ss(1) + 1.d0)
      nr(3) =  1.0d0 - ss(1) *  ss(1)

      ns(1) =  0.5d0 * ss(2) * (ss(2) - 1.d0)
      ns(2) =  0.5d0 * ss(2) * (ss(2) + 1.d0)
      ns(3) =  1.0d0 - ss(2) *  ss(2)

      nt(1) =  0.5d0 * ss(3) * (ss(3) - 1.d0)
      nt(2) =  0.5d0 * ss(3) * (ss(3) + 1.d0)
      nt(3) =  1.0d0 - ss(3) *  ss(3)

      dr(1) =  ss(1) - 0.5d0
      dr(2) =  ss(1) + 0.5d0
      dr(3) = -ss(1) - ss(1)

      ds(1) =  ss(2) - 0.5d0
      ds(2) =  ss(2) + 0.5d0
      ds(3) = -ss(2) - ss(2)

      dt(1) =  ss(3) - 0.5d0
      dt(2) =  ss(3) + 0.5d0
      dt(3) = -ss(3) - ss(3)

c     Set local 3-d shape functions

      do i = 1,27
        shp(1,i) = dr(ir(i)) * ns(is(i)) * nt(it(i))
        shp(2,i) = nr(ir(i)) * ds(is(i)) * nt(it(i))
        shp(3,i) = nr(ir(i)) * ns(is(i)) * dt(it(i))
        shp(4,i) = nr(ir(i)) * ns(is(i)) * nt(it(i))
      end do ! i

      end
