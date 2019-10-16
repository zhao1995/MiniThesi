c$Id:$
      subroutine pfeap(xl,yl,siz,color,border)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Put FEAP logo on plots

c      Inputs:
c         xl,yl     - Location to place logo
c         siz       - Size of logo
c         color     - Color for plot
c         border    - Border type: <2 = fill; >1 = line

c      Outputs:
c         none      - Plot output to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata2.h'

      integer   i, ico1, ico2, is,color,border
      real*8    xl,yl,dixl,siz,size

      integer   ifx(11),ify(11),iex(13),iey(13),iln(2)
      integer   iax(12),iay(12),ipx(15),ipy(15),ixl(4)

      save

      data ifx / 0,10,45,43,18,16,26,24,14,10, 0/
      data ify / 0,50,50,40,40,30,30,20,20, 0, 0/

      data iex / 0,10,45,43,18,16,26,24,14,12,37,35, 0/
      data iey / 0,50,50,40,40,30,30,20,20,10,10, 0, 0/

      data iax / 0,20,30,40,30,23,18,26,28,14,10, 0/
      data iay / 0,50,50, 0, 0,33,20,20,10,10, 0, 0/

      data ipx / 0,10,30,40,37,29,13,15,24,26,28,26,18,10,0/
      data ipy / 0,50,50,40,24,15,15,25,25,28,38,40,40, 0,0/

      data ixl /40,80,120,165/

c     Save line type

      ico1   = ilno(1)
      ico2   = ilno(2)
      iln(1) = 0
      iln(2) = 1
      call plline(iln)

c     Plot FEAP letters

      size = 200.0/siz
      if(border.le.1) then
        is = 1
      else
        is = 3
      endif

c     F

      call pppcol(color,0)
      dixl = xl*size + ixl(1)
      call dplot((ifx(1)+dixl)/size,ify(1)/size+yl,is)
      do i = 2,11
        call dplot((ifx(i)+dixl)/size,ify(i)/size+yl,2)
      end do ! i
      if(is.eq.1) call clpan

c     E

      dixl = xl*size + ixl(2)
      call dplot((iex(1)+dixl)/size,iey(1)/size+yl,is)
      do i = 2,13
        call dplot((iex(i)+dixl)/size,iey(i)/size+yl,2)
      end do ! i
      if(is.eq.1) call clpan

c     A

      dixl = xl*size + ixl(3)
      call dplot((iax(1)+dixl)/size,iay(1)/size+yl,is)
      do i = 2,12
        call dplot((iax(i)+dixl)/size,iay(i)/size+yl,2)
      end do ! i

      if(is.eq.1) then
        call clpan
      endif

c     P

      call pppcol(color,0)
      dixl = xl*size + ixl(4)
      call dplot((ipx(1)+dixl)/size,ipy(1)/size+yl,is)
      do i = 2,15
        call dplot((ipx(i)+dixl)/size,ipy(i)/size+yl,2)
      end do ! i

      if(is.eq.1) then
        call clpan
      endif

c     Restore line type

      ilno(1) = ico1
      ilno(2) = ico2

      end
