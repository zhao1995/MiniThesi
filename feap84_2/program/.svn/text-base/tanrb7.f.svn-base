c$Id:$
      subroutine tanrb7(mass,pi1,theta,dt,inert1,rotv,rota,nst, srb)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Tangent for rigid body: Update method 7

c      Inputs:
c         mass     - Mass of rigid body
c         pi1      - Momentum at t_n+1
c         theta(3) - Time integration parameters: beta, gamma, alpha
c         dt       - Time increment
c         inert1   - Spatial inertia at t_n+1
c         rotv(3)  - Rotational velocity
c         rota(3)  - Rotational acceleration
c         nst      - Dimension of tangent matrix, s

c      Outputs:
c         s(nst,*) - Tangent matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nst, i, j
      real*8    s(6,6),srb(nst,*), theta(3), dt, dtr, cons, mass
      real*8    pi1(3),rotv(3),rota(3),inert1(3,3),tl(3,3),h(3,3)
      real*8    temp(3,3), term(3,3)

      save

      call tmat7(rotv,rota, tl,h)

      cons = theta(2)/(theta(1)*dt)

      do i = 1,3
        do j = 1,3
          inert1(i,j) = cons*inert1(i,j)
        end do ! j

        temp(i,1) = inert1(i,2)*rotv(3) - inert1(i,3)*rotv(2)
        temp(i,2) = inert1(i,3)*rotv(1) - inert1(i,1)*rotv(3)
        temp(i,3) = inert1(i,1)*rotv(2) - inert1(i,2)*rotv(1)

      end do ! i

      temp(1,2) = temp(1,2) + pi1(3)
      temp(2,1) = temp(2,1) - pi1(3)
      temp(1,3) = temp(1,3) - pi1(2)
      temp(3,1) = temp(3,1) + pi1(2)
      temp(2,3) = temp(2,3) + pi1(1)
      temp(3,2) = temp(3,2) - pi1(1)

      do i = 1,3
        do j = 1,3
          term(i,j) = temp(i,1)*h(1,j) + temp(i,2)*h(2,j)
     &              + temp(i,3)*h(3,j) + inert1(i,j)
        end do ! j
      end do ! i

      dtr  = 1.d0/(theta(3)*dt)
      cons = mass*cons*dtr

      do i = 1,3
        do j = 1,3
          s(i  ,j  ) = 0.0d0
          s(i  ,j+3) = 0.0d0
          s(i+3,j  ) = 0.0d0
        end do ! j

        s(i,i)=s(i,i) + cons

        do j = 1,3
          s(i+3,j+3) = dtr*(term(i,1)*tl(1,j)
     &                    + term(i,2)*tl(2,j)
     &                    + term(i,3)*tl(3,j))
        end do ! j
      end do ! i

c     2-D case

      if(nst.eq.3) then
        do j = 1,2
          do i = 1,2
            srb(i,j) = s(i,j)
          end do ! i
          srb(3,j) = s(6,j)
          srb(j,3) = s(j,6)
        end do ! j
        srb(3,3) = s(6,6)

c     3-D case

      else
        do j = 1,6
          do i = 1,6
            srb(i,j) = s(i,j)
          end do ! i
        end do ! j
      endif

      end
