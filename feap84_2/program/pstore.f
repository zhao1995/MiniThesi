c$Id:$
      subroutine pstore(ni,nv,ixl,vall,ixm,valm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Move integer and/or real (double precision) array
c               to new locations.

c      Inputs:
c         ni        - Length of integer array
c         nv        - Length of real array
c         ixl(*)    - Integer array to move
c         vall(*)   - Real array to move

c      Outputs:
c         ixm(*)    - Moved integer array
c         valm(*)   - Moved real array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ni,nv, i
      integer   ixl(ni),ixm(ni)
      real*8    vall(nv),valm(nv)

      save

      do i = 1,ni
        ixm(i) = ixl(i)
      end do ! i

      do i = 1,nv
        valm(i) = vall(i)
      end do ! i

      end
