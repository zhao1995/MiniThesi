c$Id:$
      subroutine quanrm ( qua )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Normalize Unit Quaternion (vector,scalar):

c      Inputs:
c         qua(4)  - Quaternion to normalize

c      Outputs:
c         qua(4)  - Normalized quaternion
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i
      real*8    qtnrm, qua(4)

      save

      qtnrm = 1.d0/sqrt(qua(1)*qua(1) + qua(2)*qua(2)
     &                + qua(3)*qua(3) + qua(4)*qua(4))
      do i = 1,4
        qua(i) = qua(i) * qtnrm
      end do ! i

      end
