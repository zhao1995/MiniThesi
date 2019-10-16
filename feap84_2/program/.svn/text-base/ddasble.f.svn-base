c$Id:$
      subroutine dassble(s,ld,ad,ns)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Assemble diagonal coefficient array

c      Inputs:
c         s(ns,*)  - Coefficient matrix
c         ld(ns)   - Assembly locations
c         ns       - Size of matrix

c      Outputs:
c         ad(*)    - Digonal of coefficient matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    ns, ld(*)
      real*8     s(ns,ns), ad(*)

      integer    i,j

c     Loop through rows to perform assembly

      do i = 1,ns
        if(ld(i).gt.0) then

c         Loop through columns to perform assembly

          do j = 1,ns
            if(ld(j).eq.ld(i)) then
              ad(ld(i)) = ad(ld(i)) + s(i,j)
            endif
          end do ! j
        endif
      end do ! i

      end
