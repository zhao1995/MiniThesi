c$Id:$
      subroutine pblenda1(i,i1,i2,is,iside,isd,flsd)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Find side number for 2-d blend.

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

      logical   flsd
      integer   i,i1,i2,i3,i4, j,k, isd
      integer   is(isd,*),iside(*)

      save

c     Check each side

      do j = 1,numsd
        i3 = is(1,j)
        if(i3.eq.0 .or. i3.eq.1 .or.i3.eq.3) then
          k = 3
        elseif(i3.eq.2) then
          do i4 = 3,isd
            if(is(i4,j).ne.0) then
              k = i4
            endif
          end do ! i4
        endif
        if((i1.eq.is(2,j) .and. i2.eq.is(k,j)) .or.
     &     (i1.eq.is(k,j) .and. i2.eq.is(2,j))) then
           iside(i) = j
           flsd     = .false.
           return
        endif
      end do ! j
      flsd = .true.

      end
