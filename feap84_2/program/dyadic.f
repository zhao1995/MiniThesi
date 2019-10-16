c$Id:$
      subroutine dyadic(rot,norm,dyad)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute a dyadic: <dyad = eoe> and its   norm: <|O|>

c      Inputs:
c         rot   - Rotation for construction


c      Outputs:
c         norm  - Norm of rotation
c         dyad  - Dyadic for rotation
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, j
      real*8    rot(3), dyad(3,3), norm, rnorm, rotj

      save

      norm = rot(1)**2 + rot(2)**2 + rot(3)**2

      if(norm.gt.0.0d0) then

        rnorm = 1.d0/norm
        norm  = sqrt(norm)

        do j = 1,3
          rotj = rot(j)*rnorm
          do i = 1,3
            dyad(i,j) = rot(i)*rotj
          end do ! i
        end do ! j

      else

        do j = 1,3
          do i = 1,3
            dyad(i,j) = 0.0d0
          end do ! i
        end do ! j

      end if

      end
