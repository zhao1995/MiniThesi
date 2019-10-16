c$Id:$
      subroutine rotqua ( rot , qua )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Subroutine to convert 3x1 rotation vector rot to
c               4x1 unit quaternion. Added by JCS
c               Quaternions are Stored as: (vector,scalar).

c      Inputs:
c         rot(3)   - Rotation vector

c      Outputs:
c         qua(4)   - Quaternion
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i
      real*8    rotnr2, fac2, rot(3) , qua(4)

      save

      rotnr2 = rot(1)*rot(1)+rot(2)*rot(2)+rot(3)*rot(3)
      rotnr2 = sqrt(rotnr2)*0.5d0
      if (rotnr2.lt.1.d-04) then
        fac2   = 0.5d0 - rotnr2**2*(840.d0
     &                 - rotnr2**2*(42.d0
     &                 - rotnr2**2))/10080.d0
      else
        fac2   = sin(rotnr2)/rotnr2 * 0.5d0
      endif

      qua(4)    = cos(rotnr2)
      do i  = 1 , 3
        qua(i) = fac2 * rot(i)
      end do ! i

      end
