c$Id:$
      subroutine jplane(lam,rlam,ebig,alp1,x1,x2h,x2,
     &                  r1h,r1,r2h,r2,n1,pr,sr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form Planar joint arrays for connecting rigid bodies.

c      Inputs:
c         lam(3)    - Value of lagrange multiplier force on joint
c         rlam(3,-) - Array of lambda values
c         ebig(3,3) - Orientation of plane in reference frame
c         alp1      - alpha value for t_n+alpha
c         x1(3)     - Location of body 1 connection point
c         x2h(3)    - Location of body 2 connection point
c         x2(3)     - Location of body 2 connection point
c         r1(3)     - Distance from body 1 mass center to x1 at t_n+1
c         r1h(3)    - Distance from body 1 mass center to x1 at t_n+1
c         r2(3)     - Distance from body 2 mass center to x2 at t_n+1
c         r2h(3)    - Distance from body 2 mass center to x2 at t_n+1
c         n1        - Rigid body number for plane

c      Outputs:
c         pr(*)     - Residual array for joint
c         sr(*,*)   - Tangent  array for joint
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,k, n1
      real*8    alp1, alpn, lam(1), ebig(3,3), rlam(9,6,*)
      real*8    p(3,6), s(17,17), pr(3,7), sr(21,21)
      real*8    esmn(3), esm1(3), esma(3)
      real*8    x1(3),x2h(3),x2(3),r1h(3),r1(3),r2h(3),r2(3)

      save

c     Form updated basis vectors at t_n

      call quavec(rlam(1,1,n1),ebig(1,3),esmn)

c     Form updated basis vectors at t_n+1

      call quavec(rlam(1,3,n1),ebig(1,3),esm1)

c     Interpolate the base vectors to t_n+a
      alpn = 1.d0 - alp1
      do i = 1,3
        esma(i) = alpn*esmn(i) + alp1*esm1(i)
      end do ! i

      call jtmix(esm1,esma,x1,x2,r1,r2,x2h,r1h,r2h,lam(1),p,s)

      do j = 1,13
        do k = 1,13
          sr(j,k) = sr(j,k) + s(j,k)
        end do ! k
      end do ! j
      do j = 1,3
        do k = 1,4
          pr(j,k) = pr(j,k) + p(j,k)
        end do ! k
      end do ! j
      pr(1,5) = pr(1,5) + p(1,5)

      end
