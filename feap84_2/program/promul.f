c$Id:$
      subroutine promul(al,au,ad,b,c,jp,ne,add)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Routine to form c = c +/- a*b where a is square matrix
c               stored in profile form, b,c are vectors, and jp locates
c               the bottom of columns or rows in a.

c               If add .true.  : add to c,
c               If add .false. : subtract from c

c      Inputs:
c         al(*)  - Lower part of matrix
c         al(*)  - Upper part of matrix
c         ad(*)  - Diagonal of matrix
c         b(*)   - Vector to multiply
c         jp(*)  - Pointer for row/column ends of profile
c         ne     - Number equations
c         add    - Flag, Add matrix if true; else subtract

c      Outputs:
c         c(*)   - Vector with added matrix vector sum
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   add
      integer   ne, neq, jd, js, jh, j, jp(*)
      real*8    dot,ab,bj, al(*),au(*),ad(*),b(*),c(*)

      save

c     Initialize

      neq = abs(ne)
      jd  = 0

c     Multiply lower triangular part

      if(ne.gt.0) then
        do j = 1,neq
          js = jd
          jd = jp(j)
          if(jd.gt.js) then
            jh = jd - js
            ab = dot(al(js+1),b(j-jh),jh)
            if(add) then
              c(j) = c(j) + ab
            else
              c(j) = c(j) - ab
            endif
          endif
        end do ! j

c       Do diagonal part

        do j = 1,neq
          if(add) then
            c(j) = c(j) + ad(j)*b(j)
          else
            c(j) = c(j) - ad(j)*b(j)
          endif
        end do ! j
      endif

c     Multiply upper triangular part

      do j = 1,neq
        js = jd
        jd = jp(j)
        bj = b(j)
        if(add) bj = -bj
        if(jd.gt.js) then
           jh = jd - js
           call colred(au(js+1),bj,jh, c(j-jh))
        endif
      end do ! j

      end
