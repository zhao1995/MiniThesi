c$Id:$
      subroutine surf9( xi, shp, dshp, d2shp )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor           March 25, 1997            1.0

c      Acronym: Surface with 9 nodes for blend
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8   x1p,x1m,x1s,x2p,x2m,x2s,xi(2),shp(9),dshp(2,9),d2shp(3,9)

      x1p = 0.5d0*(1.d0 + xi(1))
      x1m = 0.5d0*(1.d0 - xi(1))
      x1s = (1.d0 - xi(1)*xi(1))

      x2p = 0.5d0*(1.d0 + xi(2))
      x2m = 0.5d0*(1.d0 - xi(2))
      x2s = (1.d0 - xi(2)*xi(2))

      shp(1) = x1m*x2m
      shp(2) = x1p*x2m
      shp(3) = x1p*x2p
      shp(4) = x1m*x2p

      shp(5) = x1s*x2m
      shp(6) = x1p*x2s
      shp(7) = x1s*x2p
      shp(8) = x1m*x2s

      shp(9) = x1s*x2s

      dshp(1,1) = -0.5d0*x2m
      dshp(1,2) =  0.5d0*x2m
      dshp(1,3) =  0.5d0*x2p
      dshp(1,4) = -0.5d0*x2p

      dshp(1,5) = -2.0d0*x2m*xi(1)
      dshp(1,6) =  0.5d0*x2s
      dshp(1,7) = -2.0d0*x2p*xi(1)
      dshp(1,8) = -0.5d0*x2s

      dshp(1,9) = -2.0d0*x2s*xi(1)

      dshp(2,1) = -0.5d0*x1m
      dshp(2,2) = -0.5d0*x1p
      dshp(2,3) =  0.5d0*x1p
      dshp(2,4) =  0.5d0*x1m

      dshp(2,5) = -0.5d0*x1s
      dshp(2,6) = -2.0d0*x1p*xi(2)
      dshp(2,7) =  0.5d0*x1s
      dshp(2,8) = -2.0d0*x1m*xi(2)

      dshp(2,9) = -2.0d0*x1s*xi(2)

      d2shp(1,5) = -2.d0*x2m
      d2shp(1,7) = -2.d0*x2p
      d2shp(1,9) = -2.d0*x2s

      d2shp(2,6) = -2.d0*x1p
      d2shp(2,8) = -2.d0*x1m
      d2shp(2,9) = -2.d0*x1s

      d2shp(3,1) =  0.25d0
      d2shp(3,2) = -0.25d0
      d2shp(3,3) =  0.25d0
      d2shp(3,4) = -0.25d0

      d2shp(3,5) =  xi(1)
      d2shp(3,6) = -xi(2)
      d2shp(3,7) = -xi(1)
      d2shp(3,8) =  xi(2)

      d2shp(3,9) = 4.d0*xi(1)*xi(2)

      end
