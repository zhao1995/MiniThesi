c$Id:$
      subroutine jtrots(ea,eb,eah,ebh,lambda,p,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form rotaional constraint to connect rigid bodies.

c               Rotational constraint:  ea dot eb = 0

c      Inputs:
c         ea(3)    - Orientation of axis of body a at t_n+1
c         eb(3)    - Orientation of axis of body b at t_n+1
c         eah(3)   - Orientation of axis of body a at t_n+alpha
c         ebh(3)   - Orientation of axis of body b at t_n+alpha
c         lambda   - Value of lagrange multiplier on constraint

c      Outputs:
c         p(*)     - Residual array for constraint
c         s(*,*)   - Tangent  array for constraint
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,nl
      real*8    lambda, lambdah, edot, edota, edotb
      real*8    eah(3), ebh(3), ea(3),eb(3), p(3,6), s(17,17)
      real*8    ec(3), ech(3), tlama(3,3), tlamb(3,3)

      save

      call pzero(p,18)
      call pzero(s,289)

      lambdah = 0.5d0*lambda

c     Form basic vectors:

      ec(1) = ea(2)*eb(3) - ea(3)*eb(2)
      ec(2) = ea(3)*eb(1) - ea(1)*eb(3)
      ec(3) = ea(1)*eb(2) - ea(2)*eb(1)

      ech(1) = eah(2)*ebh(3) - eah(3)*ebh(2)
      ech(2) = eah(3)*ebh(1) - eah(1)*ebh(3)
      ech(3) = eah(1)*ebh(2) - eah(2)*ebh(1)

      edot  = ea(1)*eb(1)  + ea(2)*eb(2)  + ea(3)*eb(3)
      edota = ea(1)*ebh(1) + ea(2)*ebh(2) + ea(3)*ebh(3)
      edotb = eah(1)*eb(1) + eah(2)*eb(2) + eah(3)*eb(3)

c     Form residual

      do i = 1,3
        p(i,1) = -lambda*ech(i)
        p(i,2) =  lambda*ech(i)
      end do ! i
      p(1,3) = -edot

      do i = 1,3
        do j = 1,3
          tlama(i,j) = ea(i)*ebh(j)
          tlamb(i,j) = eb(i)*eah(j)
        end do ! j
        tlama(i,i) = (tlama(i,i) - edota)
        tlamb(i,i) = (tlamb(i,i) - edotb)
      end do ! i

      nl = 7

c     Form Tangent

      do i = 1,3
        s(i  ,nl) =  ech(i)
        s(i+3,nl) = -ech(i)
        s(nl,i  ) =  ec(i)
        s(nl,i+3) = -ec(i)

        do j = 1,3
          s(i  ,j  ) =  tlama(i,j)*lambdah
          s(i+3,j  ) = -tlama(i,j)*lambdah
          s(i  ,j+3) = -tlamb(i,j)*lambdah
          s(i+3,j+3) =  tlamb(i,j)*lambdah
        end do ! j

      end do ! i

      end
