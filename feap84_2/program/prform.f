c$Id:$
      subroutine prform(neq,dr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Output of residual vector

c     Inputs:
c        neq    - Number of equations
c        dr(*)  - Vector entries
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'iodata.h'

      integer    i,neq
      real*8     dr(*)

      do i = 1,neq
        if(dr(i).ne.0.0d0) then
          write(ios,2001) i, 1, dr(i)
        endif
      end do !

c     Output last entry if zero (helps Matlab size array)

      if(dr(neq).eq.0.0d0) then
        write(ios,2001) neq, 1, dr(neq)
      endif

c     format

2001  format(2i10,1p,1d25.15)

      end
