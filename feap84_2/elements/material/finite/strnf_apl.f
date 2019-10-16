c$Id:$
      subroutine strnf_apl(f,sf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    12/12/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Tranformation Green_Cauchy strain rate to spatial
c               deformation rate

c      Inputs:
c        f(3,3)    - Deformation gradient

c      Outputs:
c        sf(3,3)   - Spatial deformation rate
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      real*8     f(3,3),sf(6,6)

      sf(1,1) = f(1,1)*f(1,1)
      sf(1,2) = f(1,2)*f(1,2)
      sf(1,3) = f(1,3)*f(1,3)
      sf(1,4) = f(1,1)*f(1,2) + f(1,2)*f(1,1)
      sf(1,5) = f(1,2)*f(1,3) + f(1,3)*f(1,2)
      sf(1,6) = f(1,1)*f(1,3) + f(1,3)*f(1,1)

      sf(2,1) = f(2,1)*f(2,1)
      sf(2,2) = f(2,2)*f(2,2)
      sf(2,3) = f(2,3)*f(2,3)
      sf(2,4) = f(2,1)*f(2,2) + f(2,2)*f(2,1)
      sf(2,5) = f(2,2)*f(2,3) + f(2,3)*f(2,2)
      sf(2,6) = f(2,1)*f(2,3) + f(2,3)*f(2,1)

      sf(3,1) = f(3,1)*f(3,1)
      sf(3,2) = f(3,2)*f(3,2)
      sf(3,3) = f(3,3)*f(3,3)
      sf(3,4) = f(3,1)*f(3,2) + f(3,2)*f(3,1)
      sf(3,5) = f(3,2)*f(3,3) + f(3,3)*f(3,2)
      sf(3,6) = f(3,1)*f(3,3) + f(3,3)*f(3,1)

      sf(4,1) = f(1,1)*f(2,1)
      sf(4,2) = f(1,2)*f(2,2)
      sf(4,3) = f(1,3)*f(2,3)
      sf(4,4) = f(1,1)*f(2,2) + f(1,2)*f(2,1)
      sf(4,5) = f(1,2)*f(2,3) + f(1,3)*f(2,2)
      sf(4,6) = f(1,1)*f(2,3) + f(1,3)*f(2,1)

      sf(5,1) = f(2,1)*f(3,1)
      sf(5,2) = f(2,2)*f(3,2)
      sf(5,3) = f(2,3)*f(3,3)
      sf(5,4) = f(2,1)*f(3,2) + f(2,2)*f(3,1)
      sf(5,5) = f(2,2)*f(3,3) + f(2,3)*f(3,2)
      sf(5,6) = f(2,1)*f(3,3) + f(2,3)*f(3,1)

      sf(6,1) = f(1,1)*f(3,1)
      sf(6,2) = f(1,2)*f(3,2)
      sf(6,3) = f(1,3)*f(3,3)
      sf(6,4) = f(1,1)*f(3,2) + f(1,2)*f(3,1)
      sf(6,5) = f(1,2)*f(3,3) + f(1,3)*f(3,2)
      sf(6,6) = f(1,1)*f(3,3) + f(1,3)*f(3,1)

      end
