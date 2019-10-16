c$Id:$
      subroutine shap1d( xi, nel, shp )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: One dimensional shape functions & natural derivatives

c      Inputs:
c         xi        - Isoparametric coordinate: ( -1 < xi < 1 )
c         nel       - Number of nodes / element   : ( 2 to 4 )

c      Outputs:
c         shp(2,3)  - Shape functions and spatial derivatives
c                     (natural derivatives only)
c         shp(1,i)  - Shape function spatial derivative: N_i,xi
c                     (natural derivatives only)
c         shp(2,i)  - Shape function                   : N_i
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nel
      real*8    xi, xi2, xi3, shp(2,nel)

      save

c     2-node shape functions and derivatives

      if(nel.eq.2) then

        shp(1,1) = -0.5d0
        shp(1,2) =  0.5d0

        shp(2,1) =  0.5d0 - 0.5d0*xi
        shp(2,2) =  0.5d0 + 0.5d0*xi

c     3-node shape functions and derivatives

      elseif(nel.eq.3) then

        xi2      =  xi*xi

        shp(1,1) =  xi - 0.5d0
        shp(1,2) =  xi + 0.5d0
        shp(1,3) = -xi - xi

        shp(2,1) =  0.5d0*(xi2 - xi)
        shp(2,2) =  0.5d0*(xi2 + xi)
        shp(2,3) =  1.0d0 - xi2

      elseif(nel.eq.4) then

        xi2      = xi*xi
        xi3      = xi*xi2

        shp(1,1) = (  1.d0 + 18.d0*xi - 27.d0*xi2)*0.0625d0
        shp(1,2) = (- 1.d0 + 18.d0*xi + 27.d0*xi2)*0.0625d0
        shp(1,3) = (-27.d0 - 18.d0*xi + 81.d0*xi2)*0.0625d0
        shp(1,4) = ( 27.d0 - 18.d0*xi - 81.d0*xi2)*0.0625d0

        shp(2,1) = (- 1.d0 +       xi + 9.d0*xi2 -  9.d0*xi3)*0.0625d0
        shp(2,2) = (- 1.d0 -       xi + 9.d0*xi2 +  9.d0*xi3)*0.0625d0
        shp(2,3) = (  9.d0 - 27.d0*xi - 9.d0*xi2 + 27.d0*xi3)*0.0625d0
        shp(2,4) = (  9.d0 + 27.d0*xi - 9.d0*xi2 - 27.d0*xi3)*0.0625d0

      endif

      end
