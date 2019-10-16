c$Id:$
      subroutine shp1dt(s,shp,nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct shp(3,3) to shp(3,4)                     02/03/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute shape functions and natural derivatives
c              at natural coordinate s.
c              Linear (2 node), quadratic (3 node) or cubic (4 node)
c              element.

c     Inputs:
c       s         : natural coordinate
c       nel       : number of nodes of element

c     Outputs:
c       shp(2,nel): shape functions and derivatives at s
c                   shp(1,1 to nel): first  derivative shape functions
c                   shp(2,1 to nel): shape functions
c                   shp(3,1 to nel): second derivative shape functions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nel
      real*8    s,shp(3,nel)

c     Linear element

      if(nel.eq.2) then

        shp(1,1) = -0.5d0
        shp(1,2) =  0.5d0

        shp(2,1) =  0.5d0 - 0.5d0*s
        shp(2,2) =  0.5d0 + 0.5d0*s

        shp(3,1) =  0.0d0
        shp(3,2) =  0.0d0

c     Quadratic element

      elseif(nel.eq.3) then

        shp(1,1) =  s - 0.5d0
        shp(1,2) =  s + 0.5d0
        shp(1,3) = -2.d0*s

        shp(2,1) =  s*(s - 1.d0)*0.5d0
        shp(2,2) =  s*(s + 1.d0)*0.5d0
        shp(2,3) =  1.d0 - s*s

        shp(3,1) =  1.0d0
        shp(3,2) =  1.0d0
        shp(3,3) = -2.0d0

c     Cubic element

      elseif(nel.eq.4) then

        shp(1,1) = 0.0625d0*( 1.d0 + 18.d0*s - 27.d0*s*s)
        shp(1,2) = 0.0625d0*(-1.d0 + 18.d0*s + 27.d0*s*s)
        shp(1,3) = 0.5625d0*(-3.d0 - 2.d0*s + 9.d0*s*s)
        shp(1,4) = 0.5625d0*( 3.d0 - 2.d0*s - 9.d0*s*s)

        shp(2,1) = 0.0625d0*(9.d0*s*s - 1.d0)*(1.d0 - 3.d0*s)
        shp(2,2) = 0.0625d0*(9.d0*s*s - 1.d0)*(1.d0 + 3.d0*s)
        shp(2,3) = 0.5625d0*(1.d0 - s*s)*(1.d0 - 3.d0*s)
        shp(2,4) = 0.5625d0*(1.d0 - s*s)*(1.d0 + 3.d0*s)

        shp(3,1) =  1.1250d0 - 3.375d0*s
        shp(3,2) =  1.1250d0 + 3.375d0*s
        shp(3,3) = -1.125d0 + 10.125d0*s
        shp(3,4) = -1.125d0 - 10.125d0*s

      endif

      end
