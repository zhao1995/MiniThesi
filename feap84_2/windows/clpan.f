c$Id:$
      subroutine clpan()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Close panel plots

c      Inputs:
c         none

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdatap.h'
      include  'pdatps.h'
      include  'pdataq.h'
      include  'plflag.h'
      include  'psdat3.h'

      integer   status, vfarea

      save

c     Close panel for filled plots

      if(screfl) status = vfarea(npf,ixy)
      npf    = 0

c     Fill panel for PostScript

      if (hdcpy .and. ipan .ge. 3 ) then
        call fppspl(ipan,xp,yp)
      endif

c     Reinitialize panel counter

      ipan = 0

      end
