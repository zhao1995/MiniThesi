c$Id:$
      function  dotx(a,b,n)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: dot (scalar) product of two vectors for (a_i - b_i)**2

c      Inputs:
c         a(*)  - Vector 1
c         b(*)  - Vector 2
c         nn    - length of vectors

c      Outputs:
c         dot   - Scalar product
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer  i,n
      real*8   dotx, a(*),b(*)

      save

      dotx = 0.0d0
      do i = 1,n
        dotx = dotx + ( a(i) - b(i) )**2
      end do ! i

      end
