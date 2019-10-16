c$Id:$
      subroutine qutrvec ( qua, vc1, vc2 )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Multiply vector vc1 by unit quaternion qua transposed
c               into vector vc2. Quaternion stored as: (vector,scalar).

c      Inputs:
c         qua(4)  - Quaternion
c         vc1(3)  - Vector 1

c      Outputs:
c         vc2(3)  - Product of quaternion transposed and vector 1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, j
      real*8    q00, q01, q02, q03, q11, q12, q13, q22, q23, q33
      real*8    qua(4), t(3,3), vc1(3), vc2(3)

      save

c     Computation of some auxiliary variables

      q00 = qua(4)*qua(4)
      q01 = qua(4)*qua(1)
      q02 = qua(4)*qua(2)
      q03 = qua(4)*qua(3)
      q11 = qua(1)*qua(1)
      q12 = qua(1)*qua(2)
      q13 = qua(1)*qua(3)
      q22 = qua(2)*qua(2)
      q23 = qua(2)*qua(3)
      q33 = qua(3)*qua(3)

c     Computation of matrix

      t(1,1) = q00 + q11 - 0.5d0
      t(2,1) = q12 - q03
      t(3,1) = q13 + q02
      t(1,2) = q12 + q03
      t(2,2) = q00 + q22 - 0.5d0
      t(3,2) = q23 - q01
      t(1,3) = q13 - q02
      t(2,3) = q23 + q01
      t(3,3) = q00 + q33 - 0.5d0

      do i = 1,3
        do j = 1,3
          t(i,j) = 2.d0*t(i,j)
        end do ! j
      end do ! i

c     Multiply matrix * vc1:

      do i = 1,3
        vc2(i) = t(i,1)*vc1(1)
     &         + t(i,2)*vc1(2) + t(i,3)*vc1(3)
      end do ! i

      end
