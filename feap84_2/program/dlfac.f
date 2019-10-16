c$Id:$
      subroutine dlfac( au, al, jp, ilim )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Reduce a single column using compact elimination

c      Inputs:
c         au(*)   - Unreduced column in upper part of A
c         al(*)   - Unreduced row    in lower part of A
c         jp(*)   - Pointer array to end of row/columns
c         ilim    - Length of reductions

c      Outputs:
c         au(*)   - Reduced column in upper part of A
c         al(*)   - Reduced row    in lower part of A
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, ic, ilim, jp(*)
      real*8    au(*), al(*), dot

      save

      do i = 2,ilim
        ic    = min(jp(i)-jp(i-1), i-1)
        au(i) = au(i) - dot(al(jp(i)+1-ic),au(i-ic),ic)
      end do ! i

      end
