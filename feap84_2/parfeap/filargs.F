c$Id:$
      subroutine filargs(nargs,iopt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'iopt' argument and adjust use for parallel  03/04/2008
c       2. Change 'include/finclude' to 'finclude'          23/01/2009
c       3. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set up filename extender for input/output, etc.

c      Inputs:
c         iopt  -- Options: 1 = command line inputs (not allowed)
c                           2 = resets for parallel runs

c      Outputs:
c         File names returned in common /comfil/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include  "finclude/petsc.h"
#     include  "setups.h"

#     include  "comfil.h"

      logical   pcomp
      integer   nargs, iopt, ibeg, iend, ires

      save

c     N.B. May need to remove call to 'doargs' or convert for computer used.
c          Call to doargs not allowed in parallel version.

      if(iopt.eq.1) then

c       call doargs(finp,fout,fres,fsav,fplt,nargs)

c     Parallel

      else

        iend = index(finp,' ') - 1
        ibeg = max(1,iend - 4)
        ires = iend - ibeg + 1
        if(ntasks.gt.1) then
          if(.not.pcomp(finp(ibeg:iend),'_0001',ires)) then
            ibeg = iend + 1
            iend = iend + 5
            finp(ibeg:iend) = '_0001'
            fout            = finp
            fout(1:1)       = 'O'
            fres            = finp
            fres(1:1)       = 'R'
            fsav            = finp
            fsav(1:1)       = 'R'
            fplt            = finp
            fplt(1:1)       = 'P'
          endif

c       Scalar

        else
          if(pcomp(finp(ibeg:iend),'_0001',ires)) then
            finp(ibeg:iend) = '     '
            fout            = finp
            fout(1:1)       = 'O'
            fres            = finp
            fres(1:1)       = 'R'
            fsav            = finp
            fsav(1:1)       = 'R'
            fplt            = finp
            fplt(1:1)       = 'P'
          endif

        endif ! ntasks
      endif ! iopt

      end
