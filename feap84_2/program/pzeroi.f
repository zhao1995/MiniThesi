c$Id:$
      subroutine pzeroi(ii,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Zero integer array of data

c      Inputs:
c         nn     - Length of array

c      Outputs:
c         ii(*)  - Array with zero values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,nn
      integer   ii(nn)

      save

      do n = 1,nn
        ii(n) = 0
      end do ! n

      end
