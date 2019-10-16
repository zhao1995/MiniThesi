c$Id:$
      subroutine eyepik(button)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Pick "eye" location for perspective views using mouse

c      Inputs:
c         none

c      Outputs:
c         button    - Mouse key pressed (l or r)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata0.h'
      include  'ppers.h'

      logical   noerr
      integer   i
      character button*1
      real*8    x1,y1,y2, qmax

      save

c     Set control values

      qmax = 0.0d0
      do i = 1,3
        qmax  = max( qmax, abs(0.5d0*(vmax(i) - vmin(i))),
     &              abs(e(i) - 0.5d0*(vmax(i) + vmin(i))) )
      end do ! i
      qmax = qmax/0.06186d0

c     Pick point from screen

      write(*,*) ' '
      write(*,*) 'PICK 1 & 2 coordinate points, or exit'
      write(*,*) ' '
      call gin(x1,y1,noerr,button)

      if(button.ne.'l') return

      e(1) = (x1 - 1.047325d0)*qmax + tg(1)
      e(2) = (y1 - 0.882500d0)*qmax + tg(2)

      write(*,*) ' '
      write(*,*) 'PICK 3 coordinate point'
      write(*,*) ' '

      call gin(x1,y2,noerr,button)

      e(3) = (y2 - 0.882500d0)*qmax + tg(3)

      do i = 1,3
        eold(i) = e(i)
      end do ! i

      end
