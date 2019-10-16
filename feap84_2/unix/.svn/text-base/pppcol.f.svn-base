c$Id:$
      subroutine pppcol(icol,jsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. c1(1) and c2(1) for passing arguments            02/07/2014
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set color for plot lines, panels, text

c      Inputs:
c         iclr     - Color parameter
c         isw      - Switch:

c      Outputs:
c         none     - Set values returned in common blocks
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata2.h'
      include  'pdatps.h'
      include  'plflag.h'

      integer   icol, isw, jsw, ncolr
      real*4    c1(1), c2(1)

      save

      data      ncolr / 17 /

c     Set color of all quantities

      isw = abs(jsw)

c     Set Postscript color values

      if (hdcpy) then
        call fppcol(icol, 2, ncolr, isw)
      endif

c     Set screen colors: Limit to 0 through ncolr

      if(screfl) then
        c1(1) = float(max(0,min(ncolr,icol)))
        c2(1) = c1(1)
        call gdx11(8,c1,c2)
      endif

      end
