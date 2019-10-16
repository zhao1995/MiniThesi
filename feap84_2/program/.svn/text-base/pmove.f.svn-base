c$Id:$
      subroutine pmove(a,b,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Move real array a into real array b

c      Inputs:
c         a(*)      - Array to move
c         nn        - Length of array to move

c      Outputs:
c         b(*)      - Moved array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,nn
      real*8    a(nn),b(nn)

      save

      do n = 1,nn
        b(n) = a(n)
      end do ! n

      end
