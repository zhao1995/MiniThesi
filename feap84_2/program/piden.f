c$Id:$
      subroutine piden(d,ns,ne)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sets array to identity (1.0)

c      Inputs:
c         ns      - First entry to set
c         ne      - Last  entry to set

c      Outputs:

c         d(*)    - Array set to unity
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,ns,ne
      real*8    d(*)

      save

      do n = ns,ne
        d(n) = 1.0d0
      end do ! n

      end
