c$Id:$
      subroutine bernst1(nr,xs,side,is,ns,ndm, x)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Construct one dimensional Bernstein interpolation for
c              coords.

c     Inputs:
c       nr        - Number of increments on side
c       xs(3,*)   - Nodal values of interpolation function
c       side      - side number (check sign)
c       is(*)     - List of side nodes
c       ns        - Order of Lagrange polynomial for side
c       ndm       - Spatial dimension of mesh

c     Outputs:
c       x(ndm,ip) - Coordinates of points
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,i1, m,m1,m2, n1,n2,n3,nr,ns,ndm, side, is(*)
      real*8    xi, xid, xs(3,*), x(ndm,*), xl(33)

      save

c     Set order of using nodes on side

      if(side.gt.0) then
        m1 = 1
        m2 = 2
        n1 = 3
        n2 = ns
        n3 = 1
      elseif(side.lt.0) then
        m1 = 2
        m2 = 1
        n1 = ns
        n2 = 3
        n3 = -1
      endif

c     Reorder coordinates

      do i = 1,ndm
        x(i,   1)    = xs(i,is(m1))
        x(i,nr+1)    = xs(i,is(m2))
        xl(i       ) = xs(i,is(m1))
        xl(ndm*(ns-1)+i) = xs(i,is(m2))
        i1 = 0
        do m = n1,n2,n3
          i1           = i1 + 1
          xl(ndm*i1+i) = xs(i,is(m))
        end do ! m
      end do ! i

c     Loop through interior points compute Bernstein interpolation

      xid =  2.0d0/dble(nr)
      xi  = -1.0d0
      do m = 1,nr-1
        xi = xi + xid
        call bernshp(xi,ns,xl,ndm, x(1,m+1) )
      end do ! m

      end
