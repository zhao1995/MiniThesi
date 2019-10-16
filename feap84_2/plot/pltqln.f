c$Id:$
      subroutine pltqln(xl, npt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: plot a curved line using quadratic isoparametric mapping

c     Nodal sequence
c              o-----o-----o
c              1     3     2

c     Inputs:
c       xl(3,3)   -  Coordinates of line
c       npt       -  Number of segments in line

c     Output:
c       Graphical line to screen/file
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    npt, n, i
      real*8     xl(3,3), xx(3), shp(3), ss, ds

c     Plot quadratic line

      ds =  2.d0/dble(npt)
      ss = -1.d0
      call plotl(xl(1,1),xl(2,1),xl(3,1),3)
      do n = 1,npt
        ss = ss + ds
        shp(1) = 0.5d0*ss*(ss - 1.0d0)
        shp(2) = 1.0d0 - ss*ss
        shp(3) = 0.5d0*ss*(ss + 1.0d0)
        do i = 1,3
          xx(i) = shp(1)*xl(i,1) + shp(2)*xl(i,2) + shp(3)*xl(i,3)
        end do ! i
        call plotl(xx(1),xx(2),xx(3),2)
      end do ! n

      end ! pltqln
