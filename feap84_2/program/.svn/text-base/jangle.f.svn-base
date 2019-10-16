c$Id:$
      subroutine jangle(j1,j2,rlam,ebig,angle)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]

c      Purpose: Determine angle between rigid bodies

c      Inputs:
c         j1        - Number of rigid body 1
c         j2        - Number of rigid body 2
c         rlam(3,*) - Array of rotation quantities for rigid bodies
c         ebig(3,3) - Orientation for revolute axis in reference frame

c      Outputs:
c         angle(3)  - Angular offset between body frames (in rads)
c                   = angle(1) = t_n
c                   = angle(2) = t_n+1
c                   = angle(3) = t_n+1/2
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,j1,j2
      real*8    rlam(9,6,*),ebig(3,3), ea(3,2,3), eb(3,3)
      real*8    x(3),y(3),angle(3),cons1

      save

      do i = 1,3
        x(i) = 0.0d0
        y(i) = 0.0d0
      end do ! i

      cons1 = 180.d0/3.141592d0

      if(j1.gt.0) then
        call quavec(rlam(1,1,j1),ebig(1,1),ea(1,1,1))
        call quavec(rlam(1,3,j1),ebig(1,1),ea(1,1,2))
        call quavec(rlam(1,1,j1),ebig(1,2),ea(1,2,1))
        call quavec(rlam(1,3,j1),ebig(1,2),ea(1,2,2))
      elseif(j1.eq.0) then
        do i = 1,3
          do j = 1,2
            ea(i,1,j) = ebig(i,1)
            ea(i,2,j) = ebig(i,2)
          end do ! j
        end do ! i
      endif

      if(j2.gt.0) then
        call quavec(rlam(1,1,j2),ebig(1,1),eb(1,1))
        call quavec(rlam(1,3,j2),ebig(1,1),eb(1,2))
      elseif(j2.eq.0) then
        do i = 1,3
          eb(i,1) = ebig(i,1)
          eb(i,2) = ebig(i,1)
        end do ! i
      endif

      do i = 1,3
        ea(i,1,3) = 0.5d0*(ea(i,1,2) + ea(i,1,1))
        ea(i,2,3) = 0.5d0*(ea(i,2,2) + ea(i,2,1))
        eb(i,3)   = 0.5d0*(eb(i,2)   + eb(i,1))
      end do ! i

      do i = 1,3
        do j = 1,3
          x(i) = x(i) + eb(j,i)*ea(j,1,i)
          y(i) = y(i) + eb(j,i)*ea(j,2,i)
        end do ! j
      end do ! i

      do i = 1,3
        angle(i) = atan2(y(i),x(i))
      end do ! i

c      write(*,*)'JANGLE',n,':angle_n = ',angle(1),'rad =',
c     &angle(1)*cons1,'deg'
c      write(*,*)'JANGLE',n,':angle_1 = ',angle(2),'rad =',
c     &angle(2)*cons1,'deg'
c      write(*,*)'JANGLE',n,':angle_h = ',angle(3),'rad =',
c     &angle(3)*cons1,'deg'

      end
