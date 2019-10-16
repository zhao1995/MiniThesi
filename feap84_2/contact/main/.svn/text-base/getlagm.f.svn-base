c$Id:$
      subroutine getlagm(ilm,lnod,nlag, lagm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor           May   30, 1996            1.0

c      Acronym: GET LAGrange Multipliers

c      Purpose: Modify profile for active contacts
c               N.B. mr(np(224)) stores the IAD array
c                    hr(np( 26)) stores the increments to multipliers

c      Inputs :
c        ilm(*)     - Lagrange multiplier equation numbers
c        lnod       - Number multiplier nodes
c        nlag       - Number of lagrange multiplier equations/node

c      Outputs:
c        lagm(*)    - Updated values of multiplier
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   i,j,ii,nn, lnod,nlag, ilm(*), idlag
      real*8    lagm(*)

c     Add Lagrange multiplier equations

      nn = 0
      fp(1) = np(26) - 1 - nlag
      do i = 1,lnod
        if(ilm(i).gt.0) then
          ii = idlag(mr(np(224)),ilm(i))
          if(ii.gt.0) then
            fp(2) = fp(1) + ii                 ! solution addr.
            do j = 1,nlag
              nn       = nn + 1
              lagm(nn) = lagm(nn) + hr(fp(2)+j)
            end do ! j
          else
            do j = 1,nlag
              nn       = nn + 1
              lagm(nn) = 0.0d0
            end do ! j
          endif
        endif
      end do ! i

      end
