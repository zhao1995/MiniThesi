c$Id:$
      subroutine genclr(ndf, v, nty, numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Remove initial conditions on nodes merged by tie
c               command

c      Inputs:
c         ndf      - Number dof/node
c         nty(*)   - Nodal type
c         numnp    - Number of nodes

c      Outputs:
c         v(ndf,*) - Initial conditions with merged nodes removed
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndf,numnp, i,n, nty(numnp)
      real*8    v(ndf,numnp)

      save

      do n = 1,numnp
        if(nty(n) .lt. 0) then
          do i = 1,ndf
            v(i,n) = 0.0d0
          end do ! i
        endif
      end do ! n

      end
