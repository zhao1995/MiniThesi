c$Id:$
      subroutine quarot ( qua, rot )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Subroutine to convert 4x1 unit quaternion qua to
c               3x1 rotation vector.
c               Quaternions are Stored as: (vector,scalar).
c      Inputs:
c         qua(4)  - Quaternion

c      Outputs:
c         rot(3)  - Rotation vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i
      real*8    qtnrm, rotnr2, rotnrm, rotfac, qua(4), rot(3)

      save

      qtnrm  = sqrt(qua(1)**2 + qua(2)**2 + qua(3)**2)
      rotnr2 = asin(min(1.d0,qtnrm))
      rotnrm = rotnr2 * 2.d0

      if (qtnrm.gt.1.d-10) then
        rotfac = rotnrm / qtnrm
      else
        rotfac = 2.d0
      endif

      do i  = 1,3
        rot(i) = rotfac * qua(i)
      end do ! i

      end
