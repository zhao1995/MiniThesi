c$Id:$
      subroutine sasbtdc(a,b,c,x)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      Purpose: Skew matrix triple product
c          (skew(a) * skew (b)^t) dot c = (a dot b) * c - (a dot c) * b
c                -          -         -    -     -    -    -     -    -

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i
      real*8    adotb, adotc, a(3),b(3),c(3), x(3)

      save

      adotb = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
      adotc = a(1)*c(1) + a(2)*c(2) + a(3)*c(3)

      do i = 1,3
        x(i) = adotb*c(i) - adotc*b(i)
      end do ! i

      end
