c$Id:$
        subroutine fparsop(fname,ftemp,nch)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set PostScript filenames to receive plot data
c               N.B. Dummy routine in serial version

c      Inputs:
c         none

c      Outputs:
c         fname     - Final name of file
c         ftemp     - Name to temporarily hold data
c         nch       - Number of characters in ftemp filename
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nch
      character fname*(*), ftemp*17

c     Dummy routine

      end
