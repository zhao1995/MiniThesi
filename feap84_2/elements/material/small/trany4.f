c$Id:$
      subroutine trany4(r,t)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version *LM*                                01/08/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Input:
c       r(3,3) - rotation matrix
c     Output:
c       t(6,6) - transformation array
c     Important note:
c       Use this array only to rotate 4th order "yield-function-type"
c       tensors P(6,6) where the scalar sig^T*P*sig is invariant under
c       rotation. Use pushr4(t,t,P,P_rot,1.d0).
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j, i1(6),i2(6)
      real*8    r(3,3), t(6,6)

      data      i1 /1,2,3,1,2,3/
      data      i2 /1,2,3,2,3,1/

c     Form transformation array for 4th order anisotropic yield tensor

      do i = 1,3
        do j = 1,3
          t(i,j) =  r(i1(j),i1(i))*r(i2(j),i2(i))
        end do ! j
        do j = 4,6
          t(i,j) =  r(i1(j),i1(i))*r(i2(j),i2(i)) *2.d0
        end do ! j
      end do ! i

      do i = 4,6
        do j = 1,3
          t(i,j) =  r(i1(j),i1(i))*r(i2(j),i2(i))
        end do ! j
        do j = 4,6
          t(i,j) =  r(i1(j),i1(i))*r(i2(j),i2(i))
     &           +  r(i1(j),i2(i))*r(i2(j),i1(i))
        end do ! j
      end do ! i

      end
