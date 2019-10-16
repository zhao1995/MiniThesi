c$Id:$
      subroutine gin(xg,yg,noerr,butn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Graphical input from mouse

c      Inputs:
c         ix,iy     - Point to place cross-hair

c      Outputs:
c         ix,iy     - Screen coordinates when button pressed
c         noerr     - Flag, true if no error
c         butn      - Button on mouse pressed
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata2.h'
      include  'plflag.h'
      include  'psdat3.h'
      include  'wdata.h'

      character butn*1
      logical   noerr
      integer   ix, iy
      real*8    xg,yg

      save

      ix = nint(xg*dble(idx))
      iy = nint(yg*dble(idy))

      if(ix.eq.0 .or. iy.eq.0) then
        ix = 11000
        iy = 11000
      end if
      if(screfl) call dosgin(ix,iy,butn)

c     Set return coordinates

      if(butn.eq.'l' .or. butn.eq.'m' .or. butn.eq.'r') then
        noerr = .true.
        xg = dble(ix - 600)/dble(idx)
        yg = dble(iy - 100)/dble(idy)
      else
        noerr = .false.
      endif

      end
