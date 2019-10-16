c$Id:$
      subroutine rnodex(ra,r,x,rcg,rlam,ixt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set coordinates for rigid node

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i,ixt
      real*8   ra(3),r(3),x(3),rcg(3,11,*),rlam(3,3,6,*)

      save

      do i = 1,3
        r(i) = rlam(i,1,3,ixt)*(x(1) - rcg(1,1,ixt))
     &       + rlam(i,2,3,ixt)*(x(2) - rcg(2,1,ixt))
     &       + rlam(i,3,3,ixt)*(x(3) - rcg(3,1,ixt))
        ra(i) = r(i)
      end do ! i

      end
