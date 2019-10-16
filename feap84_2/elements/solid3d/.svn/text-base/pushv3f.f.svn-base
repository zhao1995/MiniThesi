c$Id:$
      subroutine pushv3f (fi,v)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute push forward of vector

c     Inputs:
c        fi(3,3) - Deformation gradient
c        v(3)    - Vector

c     Outputs:
c        v(3)    - Pushed vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer i
      real*8  fi(3,3),v(3),t(3)

c     Push-forward a 1 vector.

      do i =1,3
        t(i) = v(1)*fi(1,i) + v(2)*fi(2,i) + v(3)*fi(3,i)
      end do ! i

      do i = 1,3
        v(i) = t(i)
      end do ! i

      end
