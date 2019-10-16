c$Id:$
      subroutine storeh(h,s,nst,ni)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Saves static condensation part of element tangent
c               in history (enhanced strain or mixed type elements)

c      Inputs:
c         s(nst,*)  - Array to save in history
c         nst       - Dimension of s array
c         ni        - Dimension of h array

c      Outputs:
c         h(ni,*)   - Saved history array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nst,ni, i,j
      real*8    h(ni,nst),s(nst,nst)

      save

c     Move static condensed array to history

      do i = 1,ni
        do j = 1,nst
          h(i,j) = s(i,j)
        end do ! j
      end do ! i

      end
