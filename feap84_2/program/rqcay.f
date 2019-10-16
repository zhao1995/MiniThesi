c$Id:$

      subroutine rqcay ( rot , qua )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Subroutine to convert 3x1 rotation vector rot to
c               4x1 unit quaternion according to Cayley transform.
c               Quaternions are Stored as: (vector,scalar).

c      Inputs:
c         rot(3)    - Rotation vector

c      Outputs:
c         qua(4)    - Quaternion
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i
      real*8    rot(3) , qua(4), qfac

      save

      qua(4) = 1.d0 / sqrt(1.d0
     &       + 0.25d0*(rot(1)**2 + rot(2)**2 + rot(3)**2))

      qfac   = qua(4) * 0.5d0
      do i=1,3
        qua(i) = qfac * rot(i)
      end do ! i

      end
