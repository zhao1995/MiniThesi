c$Id:$
      subroutine jdisp0(lam,def,rlam,ebig,x1,x2h,x2,
     &                  r1h,r1,r2h,r2,j1,pr,sr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form DISPLACEMENT control arrays

c      Inputs:
c         lam       - Value of lagrange multiplier force on joint
c         rlam(3,*) - Rotation array for rigid bodies
c         ebig(3,3) - Orientation of slider axis
c         x1(3)     - Location of body 1 connection point at t_n+1
c         x2h(3)    - Location of body 2 connection point at t_n+1/2
c         x2(3)     - Location of body 2 connection point at t_n+1
c         r1h(3)    - Distance from body 1 mass center to x1 at t_n+1/2
c         r1(3)     - Distance from body 1 mass center to x1 at t_n+1
c         r2h(3)    - Distance from body 2 mass center to x2 at t_n+1/2
c         r2(3)     - Distance from body 2 mass center to x2 at t_n+1
c         j1        - Number of rigid body 1

c      Outputs:
c         pr(*)    - Residual array for joint
c         sr(*,*)  - Tangent  array for joint
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'prld1.h'

      integer   i,j,j1
      real*8    defp,dist,lam, def(2), rlam(9,6,*), ebig(3,3)
      real*8    x1(3),x2h(3),x2(3),r1h(3),r1(3),r2h(3),r2(3)
      real*8    ean(3),ea1(3),eah(3),p(3,6),s(17,17),pr(3,7),sr(21,21)

      save

      dist = sqrt((x1(1) - x2(1))*(x1(1) - x2(1))
     &          + (x1(2) - x2(2))*(x1(2) - x2(2))
     &          + (x1(3) - x2(3))*(x1(3) - x2(3)))

      defp = def(1)*prldv(int(def(2)))

c     Form updated basis vectors at t_n and at t_n+1 (e_a^3)

      call quavec(rlam(1,1,j1),ebig(1,3),ean)
      call quavec(rlam(1,3,j1),ebig(1,3),ea1)

c     Interpolate the base vectors to t_n+1/2

      do i = 1,3
        eah(i) = 0.5d0*(ean(i) + ea1(i))
      end do ! i

      call jtmix(ea1,eah,x1,x2,r1,r2,x2h,r1h,r2h,lam,p,s)

      do i = 1,12
        do j = 1,12
          sr(i,j) = sr(i,j) + s(i,j)
        end do ! j
        sr(i,13) = sr(i,13) + s(i,13)
        sr(13,i) = sr(13,i) + s(13,i)
      end do ! i

      do i = 1,3
        do j = 1,4
          pr(i,j) = pr(i,j) + p(i,j)
        end do ! j
      end do ! i
      pr(1,5) = pr(1,5) + p(1,5) + defp

      write(*,*)'JDISP0:lam = ',lam,' and dist = ',dist
      write(*,*)'JDISP0:p(1,5) = ',p(1,5),' and defp = ',defp

      end
