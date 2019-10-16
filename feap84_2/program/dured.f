c$Id:$
      function dured(al,au,ad,jh)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Reduce diagonal in unsymmetric triangular decomposition

c      Inputs:
c         al(*)  - Lower terms in row
c         au(*)  - Upper terms in column
c         ad(*)  - Reduced diagonals of previous equations
c         jh     - Length of row/column

c      Outputs:
c         dured  - reduced diagonal for current equation
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   j,jh
      real*8    dured,dot, al(jh),au(jh),ad(jh)

      save

c     Scale upper U vector by reciprocal diagonals D

      do j = 1,jh
        au(j) = au(j)*ad(j)
      end do ! j

c     Dot product of L * U

      dured = dot( al(1), au(1), jh)

c     Scale lower U vector by reciprocal diagonals D

      do j = 1,jh
        al(j) = al(j)*ad(j)
      end do ! j

      end
