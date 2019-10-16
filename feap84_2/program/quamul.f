c$Id:$
      subroutine quamul ( qu1, qu2, qu3 )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Subroutine to multiply 2 quaternions qu1*qu2 into qu3.
c               Quaternions are Stored as: (vector,scalar).

c      Inputs:
c         qu1(4)  - First  quaternion
c         qu2(4)  - Second quaternion

c      Outputs:
c         qu3(4)  - Quaternion product
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    qu12nr, qu1(4), qu2(4), qu3(4)

      save

      qu12nr = qu1(1)*qu2(1) + qu1(2)*qu2(2) + qu1(3)*qu2(3)

      qu3(1) = qu1(4)*qu2(1) + qu1(1)*qu2(4)
     &       + qu1(2)*qu2(3) - qu1(3)*qu2(2)
      qu3(2) = qu1(4)*qu2(2) + qu1(2)*qu2(4)
     &       + qu1(3)*qu2(1) - qu1(1)*qu2(3)
      qu3(3) = qu1(4)*qu2(3) + qu1(3)*qu2(4)
     &       + qu1(1)*qu2(2) - qu1(2)*qu2(1)
      qu3(4) = qu1(4)*qu2(4) - qu12nr

      end
