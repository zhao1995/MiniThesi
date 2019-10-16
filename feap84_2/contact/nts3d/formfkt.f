c$Id:$
      subroutine formfkt(xi,ni,niab,niabab)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Replace '1' by '1.0d0' for shape differences     29/10/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Anna Haraldsson             February 1998            1.0

c      Acronym:

c      Purpose: Compute shape functions & derivatives for facet point.

c      Inputs:
c         xi(2)     - Natural coordinates of facet point

c      Outputs:
c         ni(*)     - Shape functions at point
c         niab(*)   - 1-st derivative of shape functions at point
c         niabab(*) - 2-nd derivative of shape functions at point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      real*8     xi1m,xi1p,xi2m,xi2p,xi(2),ni(4),niab(4,2),niabab(4,2,2)

      xi1m          =  1.0d0 - xi(1)
      xi1p          =  1.0d0 + xi(1)
      xi2m          =  1.0d0 - xi(2)
      xi2p          =  1.0d0 + xi(2)

c     ni            =  Shape functions

      ni(1)         =  0.25d0*xi1m*xi2m
      ni(2)         =  0.25d0*xi1p*xi2m
      ni(3)         =  0.25d0*xi1p*xi2p
      ni(4)         =  0.25d0*xi1m*xi2p

c     niab(k,i)     =  Derivative of shape function with respect xi(i)

      niab(1,1)     = -0.25d0*xi2m
      niab(2,1)     = -niab(1,1)

      niab(3,1)     =  0.25d0*xi2p
      niab(4,1)     = -niab(3,1)

      niab(1,2)     = -0.25d0*xi1m
      niab(4,2)     = -niab(1,2)

      niab(2,2)     = -0.25d0*xi1p
      niab(3,2)     = -niab(2,2)

c     niabab(k,i,j) =  Derivative shape function with xi(i), xi(j)

      niabab(1,1,1) =  0.d0
      niabab(2,1,1) =  0.d0
      niabab(3,1,1) =  0.d0
      niabab(4,1,1) =  0.d0
      niabab(1,2,1) =  0.25d0
      niabab(2,2,1) = -0.25d0
      niabab(3,2,1) =  0.25d0
      niabab(4,2,1) = -0.25d0
      niabab(1,1,2) =  0.25d0
      niabab(2,1,2) = -0.25d0
      niabab(3,1,2) =  0.25d0
      niabab(4,1,2) = -0.25d0
      niabab(1,2,2) =  0.d0
      niabab(2,2,2) =  0.d0
      niabab(3,2,2) =  0.d0
      niabab(4,2,2) =  0.d0

      end
