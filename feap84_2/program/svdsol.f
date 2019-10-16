c$Id:$
      subroutine svdsol(ad,au,al,jp,b,u,v,w,t,neqj,neqr,neqt,energy)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Singular valued decomposition solution

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'debugs.h'

      integer   neqj,neqr,neqt,sing, i, j, jp(*)
      real*8    ad(*),au(*),al(*),b(*),energy,wtol
      real*8    u(neqj,neqj),v(neqj,neqj),w(neqj),t(neqj)

      save

      sing = 0

      do j = 1,neqj
        do i = 1,neqj
          u(i,j) = 0.d0
        end do ! i
      end do ! j

c     Build matrix u(m,n) from vectors au,al,ad

      do i = 1,neqj-1
        u(i,i) = ad(neqr+i)
        do j = 1,i
          u(i+1,j) = al(jp(neqr+i+1) - i + j)
          u(j,i+1) = au(jp(neqr+i+1) - i + j)
        end do ! j
      end do ! i
      u(neqj,neqj) = ad(neqt)

c     SVD decomposition

      call svdcmp(u,w,v,t, neqj)

      wtol = 0.d0
      do i = 1,neqj
        wtol = max(wtol,w(i))
      end do ! i

      wtol = wtol*1.0e-6
      do i = 1,neqj
        if(w(i).lt.wtol) then
          w(i) = 0.d0
          sing = sing + 1
        else
          w(i) = 1.d0/w(i)
        end if
      end do ! i

      if(debug.and.sing.gt.0) then
        write(*,*) ' Singular matrix with', sing,' equations.'
      end if

c     SVD solution

      call svbksb(u,w,v,t,b(neqr+1), neqj, energy)

      end
