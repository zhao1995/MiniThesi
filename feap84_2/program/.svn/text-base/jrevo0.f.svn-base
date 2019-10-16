c$Id:$
      subroutine jrevo0(lam,x,rlam,ebig,alp1,j1,j2,pr,sr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form one component of revolute joint arrays
c               for connecting rigid bodies.

c      Inputs:
c         lam       - Value of lagrange multiplier force on joint
c         rlam(3,*) - Array of rotation quantities for rigid bodies
c         ebig(3,3) - Orientation for revolute axis in reference frame
c         alp1      - alpha value for t_n+alpha
c         j1        - Number of rigid body 1
c         j2        - Number of rigid body 2

c      Outputs:
c         pr(*)     - Residual array for joint
c         sr(*,*)   - Tangent  array for joint
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,k,j1,j2
      real*8    alp1, alpn, lam, x(3,3), ebig(3,3), rlam(9,6,*)
      real*8    esmn(3,3), esm1(3,3), esma(3,3)
      real*8    p(3,6), s(17,17), pr(3,7), sr(21,21)

      save

c     Form updated basis vectors for body 'j1'

      if(j1.gt.0) then
        call quavec(rlam(1,1,j1),ebig(1,3),esmn(1,3))
        call quavec(rlam(1,3,j1),ebig(1,3),esm1(1,3))
      else
        do i = 1,3
          esmn(i,3) = ebig(i,3)
          esm1(i,3) = ebig(i,3)
        end do ! i
      endif

c     Form updated basis vectors for body 'j2'

      if(j2.gt.0) then
        call quavec(rlam(1,1,j2),ebig(1,1),esmn(1,1))
        call quavec(rlam(1,1,j2),ebig(1,2),esmn(1,2))
        call quavec(rlam(1,3,j2),ebig(1,1),esm1(1,1))
        call quavec(rlam(1,3,j2),ebig(1,2),esm1(1,2))
      else
        do i = 1,3
          esmn(i,1) = ebig(i,1)
          esm1(i,1) = ebig(i,1)
          esmn(i,2) = ebig(i,2)
          esm1(i,2) = ebig(i,2)
        end do ! i
      endif

c     Interpolate the base vectors to t_n+a

      alpn = 1.d0 - alp1
      do j = 1,3
        do i = 1,3
          esma(i,j) = alpn*esmn(i,j) + alp1*esm1(i,j)
        end do ! i
      end do ! j

      i = int(x(1,3))
      call jtrots(esm1(1,3),esm1(1,i),esma(1,3),esma(1,i),lam,p,s)

      do j = 1,6
        do k = 1,6
          sr(j,k) = sr(j,k) + s(j,k)
        end do ! k
        sr(j,7) = sr(j,7) + s(j,7)
        sr(7,j) = sr(7,j) + s(7,j)
      end do ! j

      do j = 1,3
        pr(j,1) = pr(j,1) + p(j,1)
        pr(j,2) = pr(j,2) + p(j,2)
      end do ! j
      pr(1,3) = pr(1,3) + p(1,3)

      end
