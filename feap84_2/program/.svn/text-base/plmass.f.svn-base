c$Id:$
      subroutine plmass(neq,ad)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Output of diagonal (lumped) mass

c     Inputs:
c        neq    - Number of equations
c        ad(*)  - Diagonal entries of mass
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'iodata.h'

      integer    i, neq
      real*8     ad(*)

      do i = 1,neq
        if(ad(i).ne.0.0d0) then
          write(ios,2001) i,i,ad(i)
        endif
      end do ! i

c     Output last entry if zero (helps Matlab size array)

      if(ad(neq).eq.0.0d0) then
        write(ios,2001) neq, neq, ad(neq)
      endif

c     format

2001  format(2i10,1p,1d25.15)

      end
