c$Id:$
      subroutine pblenda2(i,i1,i2,is,iside,isd)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Set new side number for 2-d blend.

c     Inputs:
c        is(isd,*) - Blending side supernode lists
c        iblend(*) - Blending functions parameters/sides
c        isd       - Dimension for sides array

c     Outputs:
c        iside(4)  - Side connection numbers for face
c        flsd      - false if no new side needed
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblend.h'

      integer   i,i1,i2,i3, isd
      integer   is(isd,*),iside(*)

      save

      do i3 = 1,isd
        is(i3,numsd) = 0
      end do ! i3
      is(2,numsd) = i1
      is(3,numsd) = i2
      iside(i)    = numsd

      end
