c$Id:$
      subroutine pmovei (ia,ib, nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Move integer array ia into integer array ib

c      Inputs:
c         ia(*)     - Array to move
c         nn        - Length of array to move

c      Outputs:
c         ib(*)     - Moved array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,nn, ia(*),ib(*)

      save

      do n = 1,nn
        ib(n) = ia(n)
      end do ! n

      end
