c$Id:$
      subroutine sbodyf(d, body)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute body force values

c     Inputs:
c       d(*)    - Material parameters

c     Outputs:
c       body(*) - Body force intensities
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'eldata.h'
      include   'prld1.h'

      integer    ii
      real*8     d(*), body(3)

c     Set body load levels

      do ii = 1,3
        if(int(d(73+ii)).gt.0) then
          body(ii) = d(10+ii) + prldv(int(d(73+ii)))*d(70+ii)
        else
          body(ii) = d(10+ii)*dm
        endif
      end do ! ii

      end
