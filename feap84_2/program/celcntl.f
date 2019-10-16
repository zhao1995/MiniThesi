c$Id:$
      subroutine celcntl(ix, ic)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute contact multiplier equations

c      Inputs:
c          ix     -  Element conectivity array.

c      Outputs:
c          ic     -  Element nodal numbers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'compac.h'

      integer    i,n, ix(ncen1,*), ic(*)

      save

      do n = 1,numcels
        do i = 0,ix(ncen1,n)-1
          ic(ix(ncen+1,n)+i) = ic(ix(ncen+1,n)+i) + 1
        end do ! i
      end do ! n

      end
