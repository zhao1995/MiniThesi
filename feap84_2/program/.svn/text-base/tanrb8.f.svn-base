c$Id:$
      subroutine tanrb8(mass,pi1,ctan,dt,inert1,rotv,nst1, srb)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute linearization of rigid-body momentum equation
c                Energy-Momentum algorithm

c      Inputs:
c         mass       - Mass of rigid body
c         pi1        - Momentum at t_n+1
c         ctan(3)    - Time integration parameters
c         dt         - Time increment
c         inert1     - Spatial inertia at t_n+1
c         rotv(3)    - Rotational velocity
c         nst        - Dimension of tangent matrix, s

c      Outputs:
c         srb(nst,*) - Tangent matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'sdata.h'

      integer   nst1, i, j
      real*8    s(6,6),ctan(3), dt, mass, const
      real*8    srb(nst1,*),pi1(3),rotv(3),roth(3),inert1(3,3),tl(3,3)

      save

c     Linearization of Cayley rotation vector

      do i = 1,3
        roth(i) = rotv(i)*ctan(1)
      end do ! i

      do i = 1,3
        do j = 1,3
          tl(i,j) = roth(i)*roth(j)
        end do ! j
        tl(i,i) = tl(i,i) + 1.d0
      end do ! i

      tl(1,2) = tl(1,2) + roth(3)
      tl(2,1) = tl(2,1) - roth(3)

      tl(1,3) = tl(1,3) - roth(2)
      tl(3,1) = tl(3,1) + roth(2)

      tl(2,3) = tl(2,3) + roth(1)
      tl(3,2) = tl(3,2) - roth(1)

c     Compute the rotational tangent

      const = ctan(3)

      do i = 1,3
        do j = 1,3
          s(i+3,j+3) = const*(inert1(i,1) * tl(j,1)
     &                      + inert1(i,2) * tl(j,2)
     &                      + inert1(i,3) * tl(j,3))
        end do ! j
      end do ! i

      s(4,5) = s(4,5) + pi1(3)/dt
      s(4,6) = s(4,6) - pi1(2)/dt

      s(5,4) = s(5,4) - pi1(3)/dt
      s(5,6) = s(5,6) + pi1(1)/dt

      s(6,4) = s(6,4) + pi1(2)/dt
      s(6,5) = s(6,5) - pi1(1)/dt

c     Compute the translational tangent

      const = mass*const

      do i = 1,3
        do j = 1,3
          s(i  ,j  ) = 0.0d0
          s(i  ,j+3) = 0.0d0
          s(i+3,j  ) = 0.0d0
        end do ! j
        s(i,i) = s(i,i) + const
      end do ! i

c     2-D case

      if(ndm.eq.2) then
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
