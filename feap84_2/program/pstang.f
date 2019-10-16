c$Id:$
      subroutine pstang(neq,ir,jc,ad)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Output of symmetric sparse stored tangent

c     Inputs:
c        neq    - Number of equations
c        ir(*)  - Row pointers
c        jc(*)  - Entries in each row
c        ad(*)  - Diagonal and upper part of array
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'iodata.h'

      integer    i1,i,j, neq,ir(*),jc(*)
      real*8     ad(*)

      i1 = 1
      do i = 1,neq
        do j = i1,ir(i)
          if(ad(j).ne.0.0d0) then
            write(ios,2001) i,jc(j),ad(j)
            if(i.ne.jc(j)) then
              write(ios,2001) jc(j),i,ad(j)
            endif
          endif
        end do ! j
        i1 = ir(i) + 1
      end do ! i

c     format

2001  format(2i10,1p,1d25.15)

      end
