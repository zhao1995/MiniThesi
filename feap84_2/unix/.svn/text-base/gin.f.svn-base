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
      include  'x11f.h'

      character butn*1
      logical   noerr
      integer   ix, iy, i
      real*8    xg,yg

      save

c     X11 graphics input mode

      if(screfl) call gdx11( 12, x11, y11)

      i    = nint(x11(1))
      butn = char(i)
      if(i.eq.32 .or. i.eq.108 .or. i.eq.109 .or. i.eq.114) then
        noerr = .true.
        ix    = nint(563.2*x11(2)/xx(2))
        iy    = nint(550.0*x11(3)/xx(3))
      else
        noerr = .false.
      endif

c     Set return coordinates

      xg = dble(ix)/dble(idx)
      yg = dble(iy)/dble(idy)

      end
