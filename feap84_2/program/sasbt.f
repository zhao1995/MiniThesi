c$Id:$
      subroutine sasbt(a,b, c)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Function: c = skew(a) * skew (b)^t

c      Inputs:
c         a(3)    - First  skew array vector
c         b(3)    - Second skew array vector

c      Outputs:
c         c(3)    - Product of skew matrics
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i,j
      real*8   adotb, a(3),b(3),c(3,3)

      save

      adotb = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

      do i = 1,3
        do j = 1,3
          c(j,i) = - b(j)*a(i)
        end do ! j
        c(i,i) = c(i,i) + adotb
      end do ! i

      end
