c$Id:$
      subroutine pcompress(ld,ls,r,s,nst,nql)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compress vectors and matrix for element array

c     Inputs:
c        ld(*)     - Uncompressed equation number array
c        r(*)      - Uncompressed vector array
c        s(nst,*)  - Uncompressed matrix array
c        nst       - Leading dimension for matrix array

c     Outputs:
c        ld(*)     - Compressed equation number array
c        r(*)      - Compressed vector array
c        s(nst,*)  - Compressed matrix array
c        nql       - Number of entries in compressed array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    m,n, nst,nql,ld(*),ls(*)

      real*8     r(*),s(nst,*)

c     Build compression list

      nql = 0
      do n = 1,nst
        if(ld(n).le.0) then
          ls(n)   = 0
        else
          nql     = nql + 1
          ls(nql) = n
        endif
      end do ! n

c     Compress array

      if(nql.ne.nst) then
        do n = 1,nql
          ld(n)    = ld(ls(n))
          r(n) = r(ls(n))
          do m = 1,nql
            s(n,m) = s(ls(n),ls(m))
          end do ! m
        end do ! n
      endif

      do n = nql+1,nst
        ld(n)= 0
      end do

      end ! pcompress
