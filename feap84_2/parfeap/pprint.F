c$Id:$
      subroutine pprint(string)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add a 'call flush' to get prompts on screen      24/07/2007
c       2. Change 'include/finclude' to 'finclude'          23/01/2009
c       3. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c       4. Add another call to flush to display new text    06/02/2009
c       3. Add call to pflush                               25/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Print prompt for interactive inputs

c      Inputs:
c         string - Character string to output

c      Outputs:
c         Writes prompt string to screen
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "finclude/petsc.h"
#     include   "pfeapb.h"
#     include   "setups.h"

      character  string*(*)
      integer    l,le

      call pflush(-1)   ! Flush all file buffers to get prompt
      le = len(string)
      do l = le,1,-1
        if(string(l:l).ne.' ') exit
      end do ! l
      if(pfeap_on) then
        if(rank.eq.0) then
          write(*,2001) string(1:l)
        endif
      else
        write(*,2001) string(1:l)
      endif
      call pflush(-1)   ! Flush all file buffers to get prompt

c     Format

2001  format(a,$)

      end
