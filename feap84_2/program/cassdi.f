c$Id:$
      subroutine cassdi(ad, s, ir, jc, ld, nneq, bycol,diagin,all)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Assemble a sparse matrix A into a compressed array
c                'ad', 'au', 'al'

c      Inputs:
c         s(*)    - Diagonal array - all there
c         ir(*)   - Location non-zero entries in A by profile col/rows
c         jc(*)   - Pointer array to find entries for equations
c         ld(*)   - Local/global array to map 's' into 'A'.
c         nneq    - Size of 's'
c         bycol   - Sparse storage scheme by column if .true.
c         diagin  - Includes diagonal in sparse assembly if .true.
c         all     - Stores all terms in row/column if .true.

c      Outputs:
c         ad(*)   - Diagonal part of A
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'cdata.h'

      logical  bycol,diagin,all
      integer  nneq, i,n, inz,inza, jc(*),ir(*), ld(*)
      real*8   ad(*), s(*)

      save

c     Compact diagonal assembly of matrix

      do i = 1,nneq
        n = ld(i)
        if( n.ge.1) then

          if(all) then

c           Assemble by columns

            if(bycol) then
              if(n.eq.1) then
                inz = inza(        1, jc(n), ir, n, n)
              else
                inz = inza(jc(n-1)+1, jc(n), ir, n, n)
              endif
              ad(inz) = ad(inz) + s(i)

c           Assemble by rows (no reordering)

            else
              if(n.eq.1) then
                inz     = inza(        1, jc(n), ir, n, n)
              else
                inz     = inza(jc(n-1)+1, jc(n), ir, n, n)
              endif
              ad(inz) = ad(inz) + s(i)

            endif

c         Assemble upper/lower parts by columns

          elseif(bycol) then

c           Assemble including diagonal

            if(diagin) then
              if(n.eq.1) then
                ad(1)   = ad(1) + s(i)
              else
                inz     = inza(jc(n-1)+1, jc(n), ir, n, n)
                ad(inz) = ad(inz) + s(i)
              endif
            endif

c         Assemble upper/lower parts by rows

          else
c           Assemble including diagonal
            if(diagin) then
              inz     = inza(jc(n+neq), jc(n+neq+1)-1, ir, n, n)
              ad(inz) = ad(inz) + s(i)
            endif
          endif

c         Assemble diagonal for .not.diagin cases

          if(.not.diagin) then
            ad(n) = ad(n) + s(i)
          endif

        endif

      end do ! i

      end
