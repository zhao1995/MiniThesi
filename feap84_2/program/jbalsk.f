c$Id:$
      subroutine jbalsk(lam,x1,x2,y1,y2,yh1,yh2,al1,pr,sr,alfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form Ball & Socket (Spherical) joint arrays for
c               connecting rigid bodies.

c      Inputs:
c         lam(3)   - Value of lagrange multiplier force on joint
c         x1(3)    - Location of body 1 connection point
c         x2(3)    - Location of body 2 connection point
c         y1(3)    - Distance from body 1 mass center to x1 at t_n+1
c         y2(3)    - Distance from body 2 mass center to x2 at t_n+1
c         yh1(3)   - Distance from body 1 mass center to x1 at t_n+1/2
c         yh2(3)   - Distance from body 2 mass center to x2 at t_n+1/2
c         al1      - alpha value for t_n+alpha
c         alfl     - Tangent unsymmetric if true

c      Outputs:
c         pr(*)    - Residual array for joint
c         sr(*,*)  - Tangent  array for joint
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   alfl
      integer   i,j
      real*8    ly1dot,ly2dot,al1,aln,kappa
      real*8    x1(3),x2(3),y1(3),y2(3),yh1(3),yh2(3),ya1(3),ya2(3)
      real*8    lam(3),tlam(3),ly1(3,3),ly2(3,3), pr(3,7), sr(21,21)

      save

      data      kappa /0.d0/

c     Form residual

      aln = 1.d0 - al1
      do i = 1,3

c       Lagrange multiplier forces

        pr(i,1) = -lam(i)
        pr(i,3) =  lam(i)
        pr(i,5) =  x2(i) - x1(i)
        ya1(i)  =  yh1(i)
        ya2(i)  =  yh2(i)

c       Regularization force: kappa*(x2 - x1)


c       Total force: Lagrange multiplier + regularization forces

        tlam(i) =  lam(i) + kappa*pr(i,5)
      end do ! i

      pr(1,2) = ya1(3)*tlam(2) - ya1(2)*tlam(3)
      pr(2,2) = ya1(1)*tlam(3) - ya1(3)*tlam(1)
      pr(3,2) = ya1(2)*tlam(1) - ya1(1)*tlam(2)

      pr(1,4) = ya2(2)*tlam(3) - ya2(3)*tlam(2)
      pr(2,4) = ya2(3)*tlam(1) - ya2(1)*tlam(3)
      pr(3,4) = ya2(1)*tlam(2) - ya2(2)*tlam(1)

c     Form product of skew matrices

      ly1dot = tlam(1)*y1(1) + tlam(2)*y1(2) + tlam(3)*y1(3)
      ly2dot = tlam(1)*y2(1) + tlam(2)*y2(2) + tlam(3)*y2(3)

      do i = 1,3
        do j = 1,3
          ly1(i,j) = y1(i)*tlam(j)
          ly2(i,j) = y2(i)*tlam(j)
        end do ! j
        ly1(i,i) = ly1(i,i) - ly1dot
        ly2(i,i) = ly2(i,i) - ly2dot
      end do ! i

c     Make tangent symmetric

      if(.not.alfl) then
        ly1(1,2) = 0.5d0*(ly1(1,2) + ly1(2,1))
        ly1(2,1) = ly1(1,2)
        ly1(1,3) = 0.5d0*(ly1(1,3) + ly1(3,1))
        ly1(3,1) = ly1(1,3)
        ly1(2,3) = 0.5d0*(ly1(2,3) + ly1(3,2))
        ly1(3,2) = ly1(2,3)

        ly2(1,2) = 0.5d0*(ly2(1,2) + ly2(2,1))
        ly2(2,1) = ly2(1,2)
        ly2(1,3) = 0.5d0*(ly2(1,3) + ly2(3,1))
        ly2(3,1) = ly2(1,3)
        ly2(2,3) = 0.5d0*(ly2(2,3) + ly2(3,2))
        ly2(3,2) = ly2(2,3)
      endif

c     Form tangent

      do i = 1,3

c       Lagrange multiplier terms

        sr(i   ,i+12) =  1.d0
        sr(i+12,i   ) =  1.d0

        sr(i+6 ,i+12) = -1.d0
        sr(i+12,i+6 ) = -1.d0

        do j = 1,3
          sr(i+3,j+3) =  ly1(i,j)*al1
          sr(i+9,j+9) = -ly2(i,j)*al1
        end do ! j

      end do ! i

      sr(4   ,14  ) = -ya1(3)
      sr(5   ,13  ) =  ya1(3)
      sr(4   ,15  ) =  ya1(2)
      sr(6   ,13  ) = -ya1(2)
      sr(5   ,15  ) = -ya1(1)
      sr(6   ,14  ) =  ya1(1)

      sr(14  ,4   ) = -y1(3)
      sr(13  ,5   ) =  y1(3)
      sr(15  ,4   ) =  y1(2)
      sr(13  ,6   ) = -y1(2)
      sr(15  ,5   ) = -y1(1)
      sr(14  ,6   ) =  y1(1)

      sr(10  ,14  ) =  ya2(3)
      sr(11  ,13  ) = -ya2(3)
      sr(10  ,15  ) = -ya2(2)
      sr(12  ,13  ) =  ya2(2)
      sr(11  ,15  ) =  ya2(1)
      sr(12  ,14  ) = -ya2(1)

      sr(14  ,10  ) =  y2(3)
      sr(13  ,11  ) = -y2(3)
      sr(15  ,10  ) = -y2(2)
      sr(13  ,12  ) =  y2(2)
      sr(15  ,11  ) =  y2(1)
      sr(14  ,12  ) = -y2(1)

      end
