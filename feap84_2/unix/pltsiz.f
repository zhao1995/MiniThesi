c$Id:$
      subroutine pltsiz(nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set size of text for graphics outputs
c               N.B. Does not work with all graphics devices

c      Inputs:
c         nn        - Size of text to plot

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'debugs.h'

      integer   nn

      save

c     Not available for X11

      if(debug) then
        write(*,*) ' *WARNING* Cannot resize text to',nn
      endif

      end
