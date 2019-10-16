c$Id:$
      subroutine baseld(nneq,mp, base, u, rhs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Multiple base excitations:

c      Inputs:
c         nneq      - Number dof on nodes
c         mp        - Number base modes
c         base(*)   - Base patterns
c         u(*)      - Solution

c      Outputs:
c         rhs(*)    - Solution
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  n,nneq,mp, base(*)
      real*8   u(*), rhs(*)

      save

c     Set the base displacement pattern for each group

      do n = 1,nneq
        if(base(n).eq.mp) then
          rhs(n) = u(n)
        else
          rhs(n) = 0.0d0
        endif
      end do ! n

      end
