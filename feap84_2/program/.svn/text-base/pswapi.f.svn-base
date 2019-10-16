c$Id:$
      subroutine pswapi(a,b,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Swap integer array a with integer array b

c      Inputs:
c         a(*)      - Array a to swap
c         b(*)      - Array b to swap
c         nn        - Length of array to move

c      Outputs:
c         a(*)      - Moved a array
c         b(*)      - Moved b array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    n,nn, swap
      integer    a(nn),b(nn)

      do n = 1,nn
        swap = a(n)
        a(n) = b(n)
        b(n) = swap
      end do ! n

      end
