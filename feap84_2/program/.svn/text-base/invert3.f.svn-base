c$Id:$
      subroutine invert3(a, deta)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    20/10/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3 x 3 Matrix inverse

c      Inputs:
c        a(3,3)  - Matrix to invert

c      Outputs:
c        a(3,3)  - Inverse of matrix 'a'
c        deta    - Determinant of matrix 'a'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    i, j
      real*8     a(3,3), ai(3,3), deta, deti

c     Compute determinant

      ai(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      ai(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      ai(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)

      deta = a(1,1)*ai(1,1) + a(2,1)*ai(1,2) + a(3,1)*ai(1,3)

c     Compute inverse

      deti    = 1.d0/deta
      ai(1,1) = ai(1,1)*deti
      ai(1,2) = ai(1,2)*deti
      ai(1,3) = ai(1,3)*deti

      ai(2,1) = (a(2,3)*a(3,1) - a(3,3)*a(2,1))*deti
      ai(2,2) = (a(3,3)*a(1,1) - a(1,3)*a(3,1))*deti
      ai(2,3) = (a(1,3)*a(2,1) - a(2,3)*a(1,1))*deti

      ai(3,1) = (a(2,1)*a(3,2) - a(3,1)*a(2,2))*deti
      ai(3,2) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))*deti
      ai(3,3) = (a(1,1)*a(2,2) - a(2,1)*a(1,2))*deti

c     Store inverse back in original array

      do j = 1,3
        do i = 1,3
          a(i,j) = ai(i,j)
        end do ! i
      end do ! j

      end  ! invert3
