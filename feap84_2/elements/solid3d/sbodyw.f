c$Id:$
      subroutine sbodyw(rho,omega,xx, bf,bt, tangent)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pglob1.h'

      logical    tangent
      integer    i,j
      real*8     rho,omega,rhom2,nnyy, xx(3),yy(3),bf(3),bt(3,3)

c     Apply body forces

      do i = 1,3
        yy(i) = xx(i) - gomex(i)
      end do ! i

      rhom2 = rho*omega*omega
      nnyy  = gomev(1)*yy(1) + gomev(2)*yy(2) + gomev(3)*yy(3)
      do i = 1,3
        bf(i) = bf(i) + rhom2*(yy(i) - nnyy*gomev(i))
      end do ! i

      if(tangent) then
        do j = 1,3
          do i = 1,3
            bt(i,j) = rhom2*gomev(i)*gomev(j)
          end do ! i
          bt(j,j) = bt(j,j) - rhom2
        end do ! j
      endif

      end
