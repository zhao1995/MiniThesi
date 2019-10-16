c$Id:$
      subroutine int3dn(l,lint,sg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form Nodal quadrature points & weights for 3 dimensions

c      Inputs:
c         l       - Number of points/direction

c      Outputs:
c         lint    - Total number of points
c         sg(4,*) - Array of points and weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

c     include  'eldata.h'
      include  'iofile.h'

      integer   i,l,lint, x2(27),y2(27),z2(27),w2(27)
      real*8    h, sg(4,*)

      save

      data      x2/-1, 1, 1,-1,-1, 1, 1,-1, 0, 1, 0,-1,
     &              0, 1, 0,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0, 0/
      data      y2/-1,-1, 1, 1,-1,-1, 1, 1,-1, 0, 1, 0,
     &             -1, 0, 1, 0,-1,-1, 1, 1, 0, 0,-1, 1, 0, 0, 0/
      data      z2/-1,-1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1,
     &              1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0/
      data      w2/ 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4,
     &              4, 4, 4, 4, 4, 4, 4, 4,16,16,16,16,16,16,64/

c     Set number of total points

      lint = l*l*l

c     2x2x2 integration

      if(l.eq.8) then
        do i = 1,8
          sg(1,i) = dble(x2(i))
          sg(2,i) = dble(y2(i))
          sg(3,i) = dble(z2(i))
          sg(4,i) = 1.d0
        end do ! i

c     3x3x3 integration

      elseif(l.eq.27) then
        h = 1.0d0
        h = h/27.d0
        do i = 1,27
          sg(1,i) = dble(x2(i))
          sg(2,i) = dble(y2(i))
          sg(3,i) = dble(z2(i))
          sg(4,i) = dble(w2(i))*h
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

2000  format(' *ERROR* INT3DN: Illegal element type =',i16)

      end
