c$Id:$
      subroutine colbac(u,s,d,jj)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Column reduction for backsubstitution step.

c      Inputs:
c         s(*)  - Unreduced column
c         u(*)  - Column of upper array already reduced
c         d     - Solution value for 'u' column
c         jj    - Length to reduce

c      Outputs:
c         s(*)  - Reduced column
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   j,jj
      real*8    d,u(*),s(*)

      save

      do j = 1,jj
        s(j) = s(j) - u(j)*s(jj+1)
      end do ! j
      s(jj) = s(jj)*d

      end
