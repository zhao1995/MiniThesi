c$Id:$
      subroutine pusht2(f,k0,kt,detf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Push forward 2nd rank tensor

c       kt(i,j) = f(i,k)*k0(k,l)*f(j,l)/detf

c     Inputs:
c         k0(3,3) - Material tangent
c         f(3,3)  - deformation gradient
c         detf    - determinant of deformation gradient
c     Outputs:
c         kt(3,3) - spatialtangent
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i, j
      real*8   detf, jrec, f(3,3), k0(3,3), kt(3,3), fj(3)

c     Reciprocal deformation gradient determinant

      jrec = 1.d0/detf

c     fs = k0 * f^T

      do j = 1,3
        do i = 1,3
          fj(i) = k0(i,1)*f(j,1) + k0(i,2)*f(j,2) + k0(i,3)*f(j,3)
        end do ! i

c       kt = f * fj / detf

        do i = 1,3
          kt(i,j) = (f(i,1)*fj(1) + f(i,2)*fj(2) + f(i,3)*fj(3))*jrec
        end do ! i
      end do ! j

      end
