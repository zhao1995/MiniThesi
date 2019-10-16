c$Id:$
      subroutine vectf_apl(a,f, rf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    12/12/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set orthogonal vectors

c      Inputs:
c        a(3,2)    - First two orthogonal vectors
c        f(3,3)    - Deformation graident to push vectors to current

c      Outputs:
c        rf(3,3)   - Three orthogonal vectors (in rotated frame)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    i
      real*8     a(3,2),f(3,3),rf(3,3),r(3,3)

c     Set first two vectors from input data

      r(:,1) = a(:,1)
      r(:,2) = a(:,2)

c     Compute third orthogonal vector

      r(1,3) = r(2,1)*r(3,2) - r(3,1)*r(2,2)
      r(2,3) = r(3,1)*r(1,2) - r(1,1)*r(3,2)
      r(3,3) = r(1,1)*r(2,2) - r(2,1)*r(1,2)

c     Push forward orthgonal basis vectors: Return rotated F

      do i = 1,3
        rf(i,:) = f(i,1)*r(1,:) + f(i,2)*r(2,:) + f(i,3)*r(3,:)
      end do ! i

      end
