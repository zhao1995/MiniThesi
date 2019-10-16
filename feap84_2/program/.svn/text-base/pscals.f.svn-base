c$Id:$
      subroutine pscals(ir,jc,ad,da,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Normalize and scale stiffness terms

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    neq, n, j, ir(*),jc(*)
      real*8     ad(*),da(*)

      save

c     Normalize diagonal

      do n = 1,neq
        do j = ir(n),ir(n+1)-1
          if(jc(j).eq.n) then
            if(ad(j).ne.0.0d0) then
              da(n) = 1.d0/sqrt(abs(ad(j)))
            else
              da(n) = 1.d0
            endif
            go to 100
          endif
        end do ! j
100     continue
      end do ! n

c     Scale AU

      do n = 1,neq
        do j = ir(n),ir(n+1)-1
          ad(j) = ad(j)*da(n)*da(jc(j))
        end do ! j
      end do ! n

      end
