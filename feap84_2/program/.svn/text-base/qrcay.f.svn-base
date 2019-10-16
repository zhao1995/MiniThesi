c$Id:$
      subroutine qrcay ( qua , rot )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Subroutine to convert 4x1 unit quaternion qua to
c               3x1 rotation vector according to Cayley transform.
c               Quaternions are Stored as: (vector,scalar).

c      Inputs:
c         qua(4)  - Quaternion

c      Outputs:
c         rot(3)  - Rotation vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i
      real*8    qua(4), rot(3), qfac

      save

      qfac = 2.d0/qua(4)
      do i=1,3
        rot(i) = qfac * qua(i)
      end do ! i

      end
