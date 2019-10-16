c$Id:$
      subroutine parmacio(clab,ct,ll,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'include/finclude' to 'finclude'          23/01/2009
c       2. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c       3. Add call to pflush                               25/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Parallel solution command input subprogram

c      Inputs:
c        ll         - Command number
c        nn         - Command number

c      Outputs:
c         clab(2)   - Command character array
c         ct(3,*)   - Values of command parameters
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "finclude/petsc.h"
#     include   "pfeapb.h"
#     include   "iofile.h"
#     include   "setups.h"

      logical    pcomp,errck,tinput
      integer    ll,nn
      real*8     ct(3,*)
      character  clab(2)*15,str*36

      save

c     Input command

      call pflush(-1)
      if(ior.lt.0) then
        if(rank.eq.0) then
          clab(1) = ' '
          do while( pcomp(clab(1),' ',1) )
            write(str,2005) rank,nn+1,ll-1
            write(*,'(/)')
            call pprint(str)
            errck = tinput(clab,2,ct(1,ll),3)
          end do ! while
        else
          errck = tinput(clab,2,ct(1,ll),3)
        endif
      else
        clab(1) = ' '
        do while( pcomp(clab(1),' ',1) )
          if(ior.lt.0) then
            write(str,2005) rank,nn+1,ll-1
            call pprint(str)
          endif
          errck = tinput(clab,2,ct(1,ll),3)
        end do ! while
      endif

c     Formats

2005  format('   Proc[',i3,']  List',i3,'  Macro',i3,'> ')

      end
