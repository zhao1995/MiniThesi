c$Id:$
      subroutine pscala(ad,al,au,jp,da,neqs,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Scale the tangent matrix to have unit diagonals
c              Scaling of profile matrix to have unit diagonals

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    neqs,neq, n, j,jj, jp(*)
      real*8     ad(*),al(*),au(*),da(*)

      save

c     Normalize diagonal

      do n = 1,neq
        if(ad(n).ne.0.0d0) then
          da(n) = 1.d0/sqrt(abs(ad(n)))
          ad(n) = sign(1.d0,ad(n))
        else
          da(n) = 1.d0
        endif
      end do ! n

c     Scale AU

      do n = 2,neq
        jj = n - jp(n) - 1
        do j = jp(n-1)+1,jp(n)
          au(j) = au(j)*da(n)*da(j+jj)
        end do ! j
      end do ! n

c     Scale AL

      do n = neqs+1,neq
        jj = n - jp(n) - 1
        do j = jp(n-1)+1,jp(n)
          al(j) = al(j)*da(n)*da(j+jj)
        end do ! j
      end do ! n

      end
