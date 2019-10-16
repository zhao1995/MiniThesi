c$Id:$
      subroutine pppcol(icol,jsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
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
      integer   coli,status,vstcol,vslcol

      save

      data      ncolr / 16 /

c     Set color of all quantities

      isw = abs(jsw)

c     Set Postscript color values

      if (hdcpy) then
        call fppcol(icol, 0, ncolr, isw)
      endif

c     Set screen colors: Limit to 0 through ncolr

      coli = max(0,icol)
      if (coli .gt. ncolr ) coli = mod(coli,ncolr)
      if(screfl) status = vstcol(coli)
      if(screfl) status = vslcol(coli)

      end
