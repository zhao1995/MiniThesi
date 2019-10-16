c$Id:$
      subroutine surf4( xi, shp, dshp, d2shp )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor             May 10, 1999            1.0

c      Acronym: Surface with 4 nodes for blend
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    x1p,x1m, x2p,x2m, xi(2), shp(*), dshp(2,*), d2shp(3,*)

      x1p = 0.5*(1.d0 + xi(1))
      x1m = 0.5*(1.d0 - xi(1))

      x2p = 0.5*(1.d0 + xi(2))
      x2m = 0.5*(1.d0 - xi(2))

      shp(1) = x1m*x2m
      shp(2) = x1p*x2m
      shp(3) = x1p*x2p
      shp(4) = x1m*x2p

      dshp(1,1) = -0.5d0*x2m
      dshp(1,2) =  0.5d0*x2m
      dshp(1,3) =  0.5d0*x2p
      dshp(1,4) = -0.5d0*x2p

      dshp(2,1) = -0.5d0*x1m
      dshp(2,2) = -0.5d0*x1p
      dshp(2,3) =  0.5d0*x1p
      dshp(2,4) =  0.5d0*x1m

      d2shp(1,1) =  0.00d0
      d2shp(1,2) =  0.00d0
      d2shp(1,3) =  0.00d0
      d2shp(1,4) =  0.00d0

      d2shp(2,1) =  0.00d0
      d2shp(2,2) =  0.00d0
      d2shp(2,3) =  0.00d0
      d2shp(2,4) =  0.00d0

      d2shp(3,1) =  0.25d0
      d2shp(3,2) = -0.25d0
      d2shp(3,3) =  0.25d0
      d2shp(3,4) = -0.25d0

      end
