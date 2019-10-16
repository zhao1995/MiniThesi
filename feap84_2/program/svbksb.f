c$Id:$
      subroutine svbksb(u,w,v,t, b, mp, energy)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Singular valued decomposition backsubstitution

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, j, mp
      real*8    energy, force, u(mp,*), w(*), v(mp,*), b(*), t(*)

      save

c     Perform solution for SVD method

      do j = 1,mp
        t(j) = 0.d0
        do i = 1,mp
          t(j) = t(j) + u(i,j)*b(i)
        end do ! i
        t(j) = t(j)*w(j)
      end do ! j

      do j = 1,mp
        force = b(j)
        b(j) = 0.d0
        do i = 1,mp
          b(j) = b(j) + v(j,i)*t(i)
        end do ! i
        energy = energy + force*b(j)
      end do ! j

      end
