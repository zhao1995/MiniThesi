c$Id:$
      subroutine pfung ( d , f , detf, sig, dd, Weng)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Compute Green strain from displacement gradient  06/11/2011
c          Send displacement gradient to tranr4, add flag
c       2. Replace 'tt' by 'f' in call to pushr2            02/04/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute material stress and modulus arrays for Fung
c               Pseudo-Exponential model.
c               W(E) = C * exp(E:A:E)

c      Inputs:
c        d(*)     - Material parameters from inpt2d.
c        f(3,3)   - Deformation gradient
c        detf     - Determinant of deformation gradient

c      Outputs:
c        sig(6)   - Cauchy stress
c        dd(6,6)  - Spatial modulus array
c        Weng     - Stored energy
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    i , j
      real*8     detf, Weng, Wbar2, d(*), f(3,3), sig(6), dd(6,6)
      real*8     ee (6), ss(6), sl(6), dm(6,6), tt(6,6)

c     Compute Green-Lagrange strain measures

      ee(1) = 0.5d0*(f(1,1)*f(1,1) + f(2,1)*f(2,1) + f(3,1)*f(3,1))
     &      + f(1,1)
      ee(2) = 0.5d0*(f(1,2)*f(1,2) + f(2,2)*f(2,2) + f(3,2)*f(3,2))
     &      + f(2,2)
      ee(3) = 0.5d0*(f(1,3)*f(1,3) + f(2,3)*f(2,3) + f(3,3)*f(3,3))
     &      + f(3,3)
      ee(4) =        f(1,1)*f(1,2) + f(2,1)*f(2,2) + f(3,1)*f(3,2)
     &      + f(1,2) + f(2,1)
      ee(5) =        f(1,2)*f(1,3) + f(2,2)*f(2,3) + f(3,2)*f(3,3)
     &      + f(2,3) + f(3,1)
      ee(6) =        f(1,3)*f(1,1) + f(2,3)*f(2,1) + f(3,3)*f(3,1)
     &      + f(3,1) + f(1,3)

c     Compute basic stress measures: S = A:E

      sl(1)  = d(21)*ee(1) + d(24)*ee(2) + d(26)*ee(3)
      sl(2)  = d(24)*ee(1) + d(22)*ee(2) + d(25)*ee(3)
      sl(3)  = d(26)*ee(1) + d(25)*ee(2) + d(23)*ee(3)
      sl(4)  = d(27)*ee(4)
      sl(5)  = d(28)*ee(5)
      sl(6)  = d(29)*ee(6)

c     Exponential factor

      Wbar2  = 2.d0*d(30)*exp(ee(1)*sl(1)+ee(2)*sl(2)+ee(3)*sl(3)
     &                      + ee(4)*sl(4)+ee(5)*sl(5)+ee(6)*sl(6))

c     Compute second Piola-Kirchhoff stress

      do j = 1,6
        ss(j) = Wbar2*sl(j)
        sl(j) = sl(j)*2.d0
      end do ! j

c     Compute Material tangent moduli

      do j = 1,6
        do i = 1,6
          dm(i,j) = sl(i)*ss(j)
        end do ! i
      end do ! j

      dm(1,1) = Wbar2*d(21) + dm(1,1)
      dm(1,2) = Wbar2*d(24) + dm(1,2)
      dm(1,3) = Wbar2*d(26) + dm(1,3)

      dm(2,1) = Wbar2*d(24) + dm(2,1)
      dm(2,2) = Wbar2*d(22) + dm(2,2)
      dm(2,3) = Wbar2*d(25) + dm(2,3)

      dm(3,1) = Wbar2*d(26) + dm(3,1)
      dm(3,2) = Wbar2*d(25) + dm(3,2)
      dm(3,3) = Wbar2*d(23) + dm(3,3)

      dm(4,4) = Wbar2*d(27) + dm(4,4)
      dm(5,5) = Wbar2*d(28) + dm(5,5)
      dm(6,6) = Wbar2*d(29) + dm(6,6)

c     Push to current configuration

      call tranr4(f,f,tt,.true.)
      call pushr4(tt,tt,dm,dd,detf)
      call pushr2(f,ss,sig,detf)

c     Return energy

      Weng = 0.5d0*Wbar2

      end
