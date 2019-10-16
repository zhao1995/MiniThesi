c$Id:$
      subroutine primul(al,au,ad,u,v,jp,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Multiply profile matrix by vector using L, D, U factors.

c      Inputs:
c         al(*)  - Lower triangular part (without Identity)
c         au(*)  - Upper triangular part (without Identity)
c         ad(*)  - Diagonal terms (reciprocal)
c         u(*)   - Current vector
c         jp(*)  - Pointer array for profile
c         neq    - Number of unknowns active

c      Outputs:
c         v(*)   - Matrix times u result
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    neq, j,js,jd,jh
      integer    jp(*)
      real*8     al(*),au(*),ad(*),u(*),v(*), uj, dot

      save

c     Identity times vector

      do j = 1,neq
        v(j) = u(j)
      end do ! j

c     Multiply upper triangular part

      jd = 0
      do j = 1,neq
        js = jd
        jd = jp(j)
        if(jd.gt.js) then
          uj = -u(j)
          jh = jd - js
          call colred(au(js+1),uj,jh,v(j-jh))
        endif
      end do ! j

c     Multiply by inverse of diagonal

      do j = 1,neq
        if(ad(j).ne.0.0d0) then
          u(j) = v(j)/ad(j)
        else
          u(j) = 0.0d0
        endif
      end do ! j

c     Multiply by lower triangular matrix

      jd = 0
      do j = 1,neq
        js = jd
        jd = jp(j)
        if(jd.gt.js) then
          jh = jd - js
          v(j) = u(j) + dot(al(js+1),u(j-jh),jh)
        else
          v(j) = u(j)
        endif
      end do ! j

      end
