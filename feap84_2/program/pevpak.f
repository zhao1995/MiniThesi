c$Id:$
      subroutine pevpak(xs,ns, x,n)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Remove unwanted blanks from character string

c      Inputs:
c         xs(*)   - Unpacked character string
c         ns      - length of unpacked string

c      Outputs:
c         x(*)    - Packed character string
c         n       - length of packed string
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character xs*(*), x*(*)
      integer   ns,n,i

      save

c     Routine to pack a character string

      n = 0
      do i = 1,ns
        if(xs(i:i).ne.' ') then
          n = n + 1
          x(n:n) = xs(i:i)
        endif
      end do ! i

      end
