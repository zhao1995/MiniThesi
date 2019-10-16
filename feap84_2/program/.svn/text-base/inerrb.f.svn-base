c$Id:$
      subroutine inerrb(jj,rlam,theta,inert1,rotv,rota,
     &                  pi1,pin,w1,wn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute spatial inertia tensor for rigid body
c               and extract velocity/acceleration quantities.

c      Inputs:
c         jj(3,3)   - Material inertia tensor
c         rlam      - Current orientation (Lambda)
c         theta     - Time integration parameter

c      Outputs:
c         inert1    - Spatial inertia tensor
c         rotv(3)   - Angular velocity
c         rota(3)   - Angular acceleration
c         pi1(3)    - Current momentume
c         pin(3)    - Last time step mementum
c         w1(3)     - Current velocity
c         wn(3)     - Last time step velocity
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, j
      real*8    rlam(3,3,6)
      real*8    lam1(3,3), lamn(3,3), jj(3,3), inert1(3,3), inertn(3,3)
      real*8    lam1i(3), lamni(3)
      real*8    pi1(3), pin(3), w1(3), wn(3), theta(3), rotv(3),rota(3)

      save

c     Reassign nodal values at t_n+1 and t_n to local variables
c     Form inertia tensor from d-parameters

      do i=1,3

        wn(i)  = rlam(i,2,4)
        w1(i)  = rlam(i,2,5)
        rotv(i)= rlam(i,1,5)
        rota(i)= rlam(i,1,5)*theta(3)
      end do ! i

c     Extract rotation matrix from quaternion

      call quamat(rlam(1,1,3),lam1)
      call quamat(rlam(1,1,1),lamn)

c     Compute spatial inertia diadic and angular momenta

      do i=1,3
        do j = 1,3
          lam1i(j) = lam1(i,1)*jj(1,j) + lam1(i,2)*jj(2,j)
     &             + lam1(i,3)*jj(3,j)
          lamni(j) = lamn(i,1)*jj(1,j) + lamn(i,2)*jj(2,j)
     &             + lamn(i,3)*jj(3,j)
        end do ! j

        do j = 1,3
          inert1(i,j) = lam1i(1)*lam1(j,1) + lam1i(2)*lam1(j,2)
     &                + lam1i(3)*lam1(j,3)
          inertn(i,j) = lamni(1)*lamn(j,1) + lamni(2)*lamn(j,2)
     &                + lamni(3)*lamn(j,3)
        end do ! j

        pi1(i) = inert1(i,1)*w1(1) + inert1(i,2)*w1(2)
     &         + inert1(i,3)*w1(3)
        pin(i) = inertn(i,1)*wn(1) + inertn(i,2)*wn(2)
     &         + inertn(i,3)*wn(3)

      end do ! i

      end
