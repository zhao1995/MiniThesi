c$Id:$
      subroutine ptieint(ix1,ix2,x, ir, tol)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Tie 4-node elements with common faces

c      Inputs:
c         ix1(*)     - Interface 1 face
c         ix2(*)     - Interface f face
c         x(ndm,*)   - Nodal coordinates of mesh
c         tol        - Gap tolerance on merge

c      Outputs:
c         ir(*)      - Merge list
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'

      integer    j, n,m, n1,n2, m1,m2, ne1, ne2
      integer    ix1(*),ix2(*),ir(*)
      real*8     x(ndm,*), tol

      save

      ne1 = 0
      ne2 = 0
      do n = 1,nen
        if(ix1(n).ne.0) then
          ne1 = n
        endif
        if(ix2(n).ne.0) then
          ne2 = n
        endif
      end do ! n

      n1 = ix1(ne1)
      do n = 1, ne1
        n2 = ix1(n)
        m1 = ix2(ne2)
        do m = 1, ne2
          m2 = ix2(m)
          do j = 1,ndm
            if(abs(x(j,m1) - x(j,n2)).gt.tol .or.
     &         abs(x(j,m2) - x(j,n1)).gt.tol) then

              go to 100
            endif
          end do ! j
          ir(n1) = ir(m2)
          ir(n2) = ir(m1)
 100      m1 = m2
        end do ! m
        n1 = n2
      end do ! n

      end
