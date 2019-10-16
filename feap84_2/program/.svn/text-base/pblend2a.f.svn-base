c$Id:$
      subroutine pblend2a(iblend,iside,isd)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Construct two dimensional interpolation using blending

c     Inputs:
c        iblend(*) - Blending functions parameters/sides
c        isd       - Dimension for sides array

c     Outputs:
c        iside(4)  - Side connection numbers for face
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblend.h'
      include  'comblk.h'
      include  'pointer.h'

      logical   setvar,palloc, flsd
      integer   isd,i,i1,i2
      integer   iblend(*),iside(*), lblend(4)

      save

c     Save values of vertex nodes in case pointer to iblend changes

      do i = 1,4
        lblend(i) = iblend(10+i)
      end do ! i

c     Check each side

      do i = 1,4
        i1 = lblend(i)
        i2 = lblend(mod(i,4)+1)
        call pblenda1(i,i1,i2,mr(np(162)),iside,isd, flsd)

        if(flsd) then
          numsd = numsd + 1
          setvar = palloc( 162,'BSIDE',numsd*isd,1)
          call pblenda2(i,i1,i2,mr(np(162)),iside,isd)
        endif

      end do ! i

      end
