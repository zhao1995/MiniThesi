c$Id:$
      subroutine poldir(xlm,n,x0,y0,th)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Computes cartesian coordinates for directors
c               from polar coordinate data

c      Inputs:
c         n      - Number of director to convert
c         x0     - Center of circular arc
c         y0     - Center of circular arc
c         th     - Degree to radian conversion constant


c      Outputs:
c         xlm(*) - Cartesian data
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n
      real*8    r, x0, y0, th, xlm(9,*)

      save

c     Set polar coordinates

      r        = xlm(7,n)
      xlm(7,n) = x0 + r*cos(xlm(8,n)*th)
      xlm(8,n) = y0 + r*sin(xlm(8,n)*th)

      end
