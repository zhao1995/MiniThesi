c$Id:$
      subroutine plota(xc,yc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Places a plot box to recieve dplot/splot outputs

c      Inputs:
c         xc,yc     - Coordinates to mark where to center plot

c      Outputs:
c         none      - Plot output to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,n
      real*8    xc,yc,dx,dy,x0,y0,xh,yh,x,y, fact

      save

c     Put a box around a plot

      if(xc.eq.0.5d0) then
        fact = 1.0d0
      else
        fact = 0.5d0
      endif
      dx = 0.02d0*fact
      dy = 0.08d0*fact
      x0 = xc - 0.40d0*fact
      y0 = yc - 0.40d0*fact
      xh = xc + 0.40d0*fact
      yh = yc + 0.40d0*fact

c     Add labels to lower axis to show line direction

      call pppcol(5,-1)

      call pltext(x0,y0-0.045455d0,-1,'A')
      call pltext(xh,y0-0.045455d0,-1,'B')

c     Put label on screen (color = 6, cyan)

      call pppcol(6,-1)
      x = 1.000d0
      y = 0.782d0
      call pltext(x,y,-21,' No.    Min       Max')
      do i = 1,2
        x  = x0
        y  = y0

c       Plot vertical leg

        do n = 1,10
          call dplot(x - dx,y     ,-3)
          call dplot(x     ,y     ,-2)
          call dplot(x     ,y + dy,-2)
          call dplot(x - dx,y + dy,-2)
          y = y + dy
        end do ! n

c       Plot horizontal leg

        y  = y0
        do n = 1,10
          call dplot(x     ,y - dx,-3)
          call dplot(x     ,y     ,-2)
          call dplot(x + dy,y     ,-2)
          call dplot(x + dy,y - dx,-2)
          x = x + dy
        end do ! n
        dx = -dx
        dy = -dy
        x0 =  xh
        y0 =  yh
      end do ! i

      end
