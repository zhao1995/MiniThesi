c$Id:$
      subroutine sumcnt(ic,neq,kp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sum counts to get storage of connections

c      Inputs:
c         ic(*)  - Pointer array for equation connection list
c         neq    - Number of total equations in system

c      Outputs:
c         kp     - Length of element connection to equation array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    i,kp,neq,ic(*)

      save

c     Set up pointers.

      do i = 2, neq
         ic(i) = ic(i) + ic(i-1)
      end do ! i

      kp = ic(neq)

      end
