c$Id:$
      subroutine pblend1a(is,iblend,iside,isd)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Construct one dimensional interpolation using blending

c     Inputs:
c        is(isd,*) - Blending side supernode lists
c        iblend(*) - Blending functions parameters/sides
c        isd       - Dimension for sides array

c     Outputs:
c        iside     - Number of side to construct
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cblend.h'
      include   'cdata.h'
      include   'iofile.h'
      include   'pointer.h'
      include   'region.h'
      include   'trdata.h'
      include   'comblk.h'

      logical    setvar, palloc
      integer    isd,iside(*), is(isd,*),iblend(*)
      integer    i, j, i1,i2,i3,i4

      save

c     Set side number to use

      i1 = iblend(11)
      i2 = iblend(12)
      do j = 1,numsd
        i3 = is(1,j)
        if(i3.eq.2) then
          do i4 = 3,isd
            if(is(i4,j).ne.0) then
              i = i4
            endif
          end do ! i4
        else
          i = 3
        endif
        if((i1.eq.is(2,j) .and. i2.eq.is(i,j)) .or.
     &     (i1.eq.is(i,j) .and. i2.eq.is(2,j))) then
          iside(1) = j
          return
        endif
      end do ! j

c     Add new side

      numsd  = numsd + 1
      setvar = palloc( 162,'BSIDE',numsd*isd,1)
      i3     = 1
      call pblenda2(i3,i1,i2,mr(np(162)),iside,isd)

      end
