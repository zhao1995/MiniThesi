c$Id:$
      subroutine dufac( au, al, at, jp, ilim, je, jq )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Reduce a single column using compact elimination

c      Inputs:
c         au(*)  - Unreduced column of upper part of A
c         al(*)  - Row of unsymmetric lower part of A
c         at(*)  - Row of   symmetric lower part of A
c         jp(*)  - Pointer array for ends of row/columns
c         ilim   - Limit of symmetric reduction
c         je     - Pointer of last symmetric equation
c         jq     - Limit of unsymmetric reduction

c      Outputs:
c         au(*)  - Reduced column of upper part of A
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, ic, i1, ilim, je, jq, j1, jp(*)
      real*8    au(*), al(*), at(*), dot

      save

c     Reduce any symmetric terms

      i1 = ilim - jq
      do i = 2,i1
        ic    = min(jp(i)-jp(i-1), i-1)
        au(i) = au(i) - dot(at(jp(i)+1-ic),au(i-ic),ic)
      end do ! i

c     Reduce any unsymmetric terms
      j1 = 1 - je
      do i = max(2,i1+1),ilim
        ic    = min(jp(i)-jp(i-1), i-1)
        au(i) = au(i) - dot(al(jp(i)+j1-ic),au(i-ic),ic)
      end do ! i

      end
