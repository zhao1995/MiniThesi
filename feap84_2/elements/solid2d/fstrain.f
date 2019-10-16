c$Id:$
      subroutine fstrain(f,finv, egreen, ealmansi)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute Green-Lagrange and Almansi strains for output

c      Inputs:
c         f(3,3)      - Deformation tensor
c         finv(3,3)   - Inverse deformation tensor

c      Outputs:
c         egreen(6)   - Green-Lagrange strain
c         ealmansi(6) - Almansi strain
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      real*8     f(3,3), finv(3,3), egreen(6),ealmansi(6)

c     Green-Lagrange strains

      egreen(1) = 0.5d0*(f(1,1)*f(1,1)+f(2,1)*f(2,1)+f(3,1)*f(3,1))
     &          - 0.5d0
      egreen(2) = 0.5d0*(f(1,2)*f(1,2)+f(2,2)*f(2,2)+f(3,2)*f(3,2))
     &          - 0.5d0
      egreen(3) = 0.5d0*(f(1,3)*f(1,3)+f(2,3)*f(2,3)+f(3,3)*f(3,3))
     &          - 0.5d0
      egreen(4) = f(1,1)*f(1,2)+f(2,1)*f(2,2)+f(3,1)*f(3,2)
      egreen(5) = f(1,2)*f(1,3)+f(2,2)*f(2,3)+f(3,2)*f(3,3)
      egreen(6) = f(1,3)*f(1,1)+f(2,3)*f(2,1)+f(3,3)*f(3,1)

c     Almansi strains

      ealmansi(1) = 0.5d0 - 0.5d0*(finv(1,1)*finv(1,1)
     &                           + finv(2,1)*finv(2,1)
     &                           + finv(3,1)*finv(3,1))
      ealmansi(2) = 0.5d0 - 0.5d0*(finv(1,2)*finv(1,2)
     &                           + finv(2,2)*finv(2,2)
     &                           + finv(3,2)*finv(3,2))
      ealmansi(3) = 0.5d0 - 0.5d0*(finv(1,3)*finv(1,3)
     &                           + finv(2,3)*finv(2,3)
     &                           + finv(3,3)*finv(3,3))
      ealmansi(4) = - finv(1,1)*finv(1,2)
     &              - finv(2,1)*finv(2,2)
     &              - finv(3,1)*finv(3,2)
      ealmansi(5) = - finv(1,2)*finv(1,3)
     &              - finv(2,2)*finv(2,3)
     &              - finv(3,2)*finv(3,3)
      ealmansi(6) = - finv(1,3)*finv(1,1)
     &              - finv(2,3)*finv(2,1)
     &              - finv(3,3)*finv(3,1)

      end
