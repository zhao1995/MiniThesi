c$Id:$
      subroutine colred(au,xj,nn, b)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Columnwise reduction for back substitution

c      Inputs:
c         au(*)   - Upper column of reduced array A
c         xj      - Solution of reduced column
c         nn      - Length to reduce
c         b(*)    - Unreduced column

c      Outputs:
c         b(*)    - Reduced column
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,nn
      real*8    xj, au(*),b(*)

      save

      do n = 1,nn
        b(n) = b(n) - au(n)*xj
      end do ! n

      end
