c$Id:$
      subroutine quamat ( qua, t)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Subroutine to convert quaternion to matrix.
c               Quaternions are Stored as: (vector,scalar).

c      Inputs:
c         qua(4)  - Quaternion

c      Outputs:
c         t(3,3)  - Rotation matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    q00, q01, q02, q03, q11, q12, q13, q22, q23, q33
      real*8    qua(4), t(3,3)

      save

c     Computation of some auxiliary variables

      q00 = qua(4)*qua(4)*2.d0 - 1.d0
      q01 = qua(4)*qua(1)*2.d0
      q02 = qua(4)*qua(2)*2.d0
      q03 = qua(4)*qua(3)*2.d0
      q11 = qua(1)*qua(1)*2.d0
      q12 = qua(1)*qua(2)*2.d0
      q13 = qua(1)*qua(3)*2.d0
      q22 = qua(2)*qua(2)*2.d0
      q23 = qua(2)*qua(3)*2.d0
      q33 = qua(3)*qua(3)*2.d0

c     Computation of matrix

      t(1,1) = q00 + q11
      t(2,1) = q12 + q03
      t(3,1) = q13 - q02
      t(1,2) = q12 - q03
      t(2,2) = q00 + q22
      t(3,2) = q23 + q01
      t(1,3) = q13 + q02
      t(2,3) = q23 - q01
      t(3,3) = q00 + q33

      end
