c$Id:$
      subroutine filargs(nargs,iopt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add iopt argument and use                        03/04/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Serial version to get command line arguments

c      Inputs:
c         iopt   -  Option: 1 = command line
c                           2 = reset for parallel solution only

c      Outputs:
c         File names returned in common /comfil/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'

      integer   nargs,iopt

      save

c     N.B. May need to remove call to 'doargs' or convert for computer
c          Call to doargs not allowed in parallel version.

      if(iopt.eq.1) then
        call doargs(finp,fout,fres,fsav,fplt,nargs)
      endif

      end
