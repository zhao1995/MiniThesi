c$Id:$
      subroutine shp3en(ss,xsj,shpi,xl,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Enhanced shape functions

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer   ndm, i, j, k
      real*8    xl(ndm,*), ss(3), shpi(3,*), txij(3), xs0(3,3), xsj

c     Set up interpolations

      do i = 1,3
        txij(i) = -2.0d0*ss(i)/xsj

        xs0(i,1) = 0.125d0*(-xl(i,1) + xl(i,2) + xl(i,3) - xl(i,4)
     &                      -xl(i,5) + xl(i,6) + xl(i,7) - xl(i,8))
        xs0(i,2) = 0.125d0*(-xl(i,1) - xl(i,2) + xl(i,3) + xl(i,4)
     &                      -xl(i,5) - xl(i,6) + xl(i,7) + xl(i,8))
        xs0(i,3) = 0.125d0*(-xl(i,1) - xl(i,2) - xl(i,3) - xl(i,4)
     &                      +xl(i,5) + xl(i,6) + xl(i,7) + xl(i,8))
      end do ! i

c     Compute assumed strain 'incompatible'  shape functions

      do i = 1,3
        j = mod(i,3) + 1
        k = mod(j,3) + 1
        shpi(1,i) = txij(i)*(xs0(2,j)*xs0(3,k) - xs0(3,j)*xs0(2,k))
        shpi(2,i) = txij(i)*(xs0(3,j)*xs0(1,k) - xs0(1,j)*xs0(3,k))
        shpi(3,i) = txij(i)*(xs0(1,j)*xs0(2,k) - xs0(2,j)*xs0(1,k))
      end do ! i

      end
