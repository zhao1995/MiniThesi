c$Id:$
      subroutine plttxt(icol,xs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Place specified text on graphics screen using mouse

c      Inputs:
c         xs(2)     - Positioning of mouse cursor

c      Outputs:
c         none      - Plot text output placed on screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pdatap.h'
      include  'pdatxt.h'

      character tx*60,button*1
      logical   noerr
      integer   icol,nn
      real*8    x1,y1, xs(2)

      save

c     Request text inputs in text mode

      tx = ' '
      if(ior.gt.0) then
        if(xs(1).eq.0.0d0 .and. xs(2).eq.0.0d0 ) return
        read(ior,1000) tx
      else
        call pprint(' Enter text > ')
        read(*,1000) tx
      endif

c     Compute length of text record

      do nn = 60,1,-1
         if(tx(nn:nn).ne.' ') go to 110
      end do ! nn
      if(ior.gt.0) then
        write(iow,2001)
      else
        write(*,2001)
      endif
      return

c     Enter graphics mode, position, and output text

110   noerr = .true.
      if(xs(1).eq.0.0d0 .and. xs(2).eq.0.0d0) then
        write(*,2004)
        call plopen
        call gin(x1,y1,noerr,button)
        call pppcol(icol,1)
      else
        x1 = 0.96d0*max(0.d0,min(1.4d0,xs(1)))
        y1 = 0.96d0*max(0.d0,min(1.0d0,xs(2)))
        call plopen
      endif
      if(noerr) then
        dtext = 0.0d0
        call pltsiz(nsizt)
        call tplot(x1,y1,tx,nn,0)
        call pltsiz(1)
      else
        write(*,2002)
      endif

c     Formats

1000  format(a)

2001  format(' *WARNING* no text input ')
2002  format(' *ERROR* PLTTXT: No text plotted')
2004  format(' Position ; Enter p or press left mouse button.')

      end
