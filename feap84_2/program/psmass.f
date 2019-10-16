c$Id:$
      subroutine psmass(neq,ir,jc,ad,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Output of consistent mass/damping array (sparse)

c     Inputs:
c        neq    - Number of equations
c        ir(*)  - Row pointers
c        jc(*)  - Entries in each row
c        ad(*)  - Diagonal and upper part of array
c        isw    - Switch: 1 = diagonal; 2 = symmetric; 3 = unsymmetric
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'iodata.h'

      integer    isw,i1,i,j, neq,ir(*),jc(*)
      real*8     ad(*)

c     Output diagonal entries

      do i = 1,neq
        if(ad(i).ne.0.0d0) then
          write(ios,2001) i,i,ad(i)
        endif
      end do ! i

c     Output the last element if it is zero (for Matlab dimensioning)

      if(ad(neq).eq.0.0d0) then
        write(ios,2001) neq,neq,ad(neq)
      endif

c     Output off-diagonal entries

      if(isw.ge.2) then
        i1 = 1
        do i = 1,neq
          do j = i1,ir(i)
            if(ad(j+neq).ne.0.0d0) then
              write(ios,2001) jc(j),i,ad(j+neq)
              if(isw.eq.2) then
                write(ios,2001) i,jc(j),ad(j+neq)
              endif
            endif
            if(isw.eq.3 .and. ad(j+neq+ir(neq)).ne.0.0d0) then
              if(i.ne.jc(j)) then
                write(ios,2001) i,jc(j),ad(j+neq+ir(neq))
              endif
            endif
          end do ! j
          i1 = ir(i) + 1
        end do ! i
      endif

c     format

2001  format(2i10,1p,1d25.15)

      end
