c$Id:$
      subroutine pick()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Use graphical input to clip window for plots

c      Inputs:
c         none

c      Outputs:
c         none      - Plot picks through mouse driver devide
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata1.h'
      include  'pdata4.h'
      include  'plclip.h'
      include  'plflag.h'

      logical   noerr
      character button*1
      real*8    x1,y1,x2,y2, xy(2,2)

      save

c     Pick a point from screen

      write(*,*) 'Pick first point'
      call gin(x1,y1,noerr,button)

      xy(1,1) = 0.5d0*(sx(1) + (x1 - s0(1))/scale)
      xy(2,1) = 0.5d0*(sx(2) + (y1 - s0(2))/scale)

      write(*,*) 'Pick second point'
      call gin(x2,y2,noerr,button)

      xy(1,2) = 0.5d0*(sx(1) + (x2 - s0(1))/scale)
      xy(2,2) = 0.5d0*(sx(2) + (y2 - s0(2))/scale)

      if(x1.ne.x2 .or. y1.ne.y2) then
        nzm1 = 0
        nzm2 = 0
        call frame(xy,2,2,1)

        wmin(1) = 0.105d0
        wmax(1) = 0.895d0
        wmin(2) = 0.105d0
        wmax(2) = 0.895d0

        fwin    = .true.
        clip    = .false.

        xc(1)   = real(wmin(1))
        xc(2)   = real(wmax(1))
        xc(3)   = real(wmax(1))
        xc(4)   = real(wmin(1))
        yc(1)   = real(wmin(2))
        yc(2)   = real(wmin(2))
        yc(3)   = real(wmax(2))
        yc(4)   = real(wmax(2))
      else
        write(*,*) ' Do again, same point picked'
      endif

      end
