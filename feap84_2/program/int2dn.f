c$Id:$
      subroutine int2dn(l,lint,sg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'one3' from 'pconstant.h' instead of 'third'   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form Nodal quadrature points & weights for 2 dimensions

c      Inputs:
c         l       - Number of points/direction

c      Outputs:
c         lint    - Total number of points
c         sg(3,*) - Array of points and weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

c     include  'eldata.h'
      include  'iofile.h'
      include  'pconstant.h'

      integer   i,l,lint, x2(9),y2(9),w2(9), x3(16),y3(16),w3(16)
      real*8    g,h, sg(3,*)

      save

      data      x2/-1, 1, 1,-1, 0, 1, 0,-1, 0/
      data      y2/-1,-1, 1, 1,-1, 0, 1, 0, 0/
      data      w2/ 1, 1, 1, 1, 4, 4, 4, 4,16/
      data      x3/-3, 3, 3,-3,-1, 1, 3, 3, 1,-1,-3,-3,-1, 1, 1,-1/
      data      y3/-3,-3, 3, 3,-3,-3,-1, 1, 3, 3, 1,-1,-1,-1, 1, 1/
      data      w3/ 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 9, 9, 9, 9/

c     Set number of total points

      lint = l

c     2x2 integration: 4-node element

      if(l.eq.4) then
        do i = 1,4
          sg(1,i) = dble(x2(i))
          sg(2,i) = dble(y2(i))
          sg(3,i) = 1.d0
        end do ! i

c     3x3 integration: 9-node element

      elseif(l.eq.9) then
        h = one9
        do i = 1,9
          sg(1,i) = dble(x2(i))
          sg(2,i) = dble(y2(i))
          sg(3,i) = dble(w2(i))*h
        end do ! i

c     4x4 integration: 16-node element

      elseif(l.eq.16) then
        g = one3
        h = 0.0625d0
        do i = 1,16
          sg(1,i) = dble(x3(i))*g
          sg(2,i) = dble(y3(i))*g
          sg(3,i) = dble(w3(i))*h
        end do ! i

c     Error

      else

        write(ilg,2000) l
        write(iow,2000) l
        if(ior.lt.0) then
          write(*,2000) l
        endif
        call plstop()

      endif

c     Format

2000  format(' *ERROR* INT2DN: Illegal element type =',i16)

      end
