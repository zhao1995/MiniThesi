c$Id:$
      subroutine pptang(neq,jp,ad, al)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Output of profile stored tangent

c     Inputs:
c        neq    - Number of equations
c        jp(*)  - Column pointers
c        ad(*)  - Diagonal and upper part of array
c        al(*)  - Lower part of array
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'iodata.h'

      integer    ii,i,j, neq,jp(*)
      real*8     ad(*), al(*)

c     Output diagonal entries

      do i = 1,neq
        if(ad(i).ne.0.0d0) then
          write(ios,2001) i,i, ad(i)
        endif
      end do ! i

c     Output last entry if zero (helps Matlab size array)

      if(ad(neq).eq.0.0d0) then
        write(ios,2001) neq, neq, ad(neq)
      endif

c     Output off-diagonal entries

      do j = 2,neq
        ii = j - jp(j) + jp(j-1)
        do i = jp(j-1)+1,jp(j)
          if(ad(neq+i).ne.0.0d0) then
            write(ios,2001) ii,j,ad(neq+i)
          endif
          if(al(i).ne.0.0d0) then
            write(ios,2001) j,ii,al(i)
          endif
          ii = ii + 1
        end do ! i
      end do ! j

c     format

2001  format(2i10,1p,1d25.15)

      end
