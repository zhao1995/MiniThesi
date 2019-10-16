c$Id:$
      subroutine pushf2(g,s,sig,detf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Push forward 2nd rank tensor

c       sigma(i,j) = f(i,k)*s(k,l)*f(j,l)/detf

c     Inputs:
c         s(6)   - material stress (2nd Piola-Kirchhoff)
c         g(3,3) - deformation gradient minus identity
c         detf   - determinant of deformation gradient
c     Outputs:
c         sig(6) - spatial stress (Cauchy)
c                       | sig(1)  sig(4)  sig(6) |
c         sigma(i,j) =  | sig(4)  sig(2)  sig(5) |
c                       | sig(6)  sig(5)  sig(3) |
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i
      real*8   detf, jrec, g(3,3), f(3,3), s(6), sig(6), fs(3,3)

c     Reciprocal deformation gradient determinant

      jrec = 1.d0/detf

c     Form deformation gradient

      f = g
      do i = 1,3
        f(i,i) = f(i,i) + 1.0d0
      end do ! i

c     fs = f^t*s

      do i = 1,3
        fs(i,1) = f(i,1)*s(1) + f(i,2)*s(4) + f(i,3)*s(6)
        fs(i,2) = f(i,1)*s(4) + f(i,2)*s(2) + f(i,3)*s(5)
        fs(i,3) = f(i,1)*s(6) + f(i,2)*s(5) + f(i,3)*s(3)
      end do ! i

c     sig = fs*f/j

      sig(1) = (fs(1,1)*f(1,1) + fs(1,2)*f(1,2) + fs(1,3)*f(1,3))*jrec
      sig(2) = (fs(2,1)*f(2,1) + fs(2,2)*f(2,2) + fs(2,3)*f(2,3))*jrec
      sig(3) = (fs(3,1)*f(3,1) + fs(3,2)*f(3,2) + fs(3,3)*f(3,3))*jrec
      sig(4) = (fs(1,1)*f(2,1) + fs(1,2)*f(2,2) + fs(1,3)*f(2,3))*jrec
      sig(5) = (fs(2,1)*f(3,1) + fs(2,2)*f(3,2) + fs(2,3)*f(3,3))*jrec
      sig(6) = (fs(3,1)*f(1,1) + fs(3,2)*f(1,2) + fs(3,3)*f(1,3))*jrec

      end
