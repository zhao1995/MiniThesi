c$Id:$
      subroutine pzerol(fl,val,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set an array to logical value = val

c      Inputs:
c         val       - Logical value: true or false
c         nn        - Length of array to set

c      Outputs:
c         fl(*)     - Array set to logical state val
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,nn
      logical   fl(nn),val

      save

      do n = 1,nn
        fl(n) = val
      end do ! n

      end
