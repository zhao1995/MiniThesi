c$Id:$
      subroutine sphdir(xlm,n,x0,y0,z0,sn9)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Computes cartesian from spherical coordinates

c      Inputs:
c         xlm(9,*)   - Spherical coordinates for point
c         n          - Number of points
c         x0,y0,z0   - Origin for spherical coordinate system

c      Outputs:
c         xlm(9,*)   - Cartesian coordinates for point
c         sn9        - Sine of meridian angle
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n
      real*8    x0,y0,z0,r,sn8,cn8,sn9,cn9, xlm(9,*)

      save

c     Set spherical coordinates

      call pdegree(xlm(8,n), sn8,cn8)
      call pdegree(xlm(9,n), sn9,cn9)
      r        = xlm(7,n)
      xlm(7,n) = x0 + r*cn8*sn9
      xlm(8,n) = y0 + r*sn8*sn9
      xlm(9,n) = z0 + r*cn9

      end
