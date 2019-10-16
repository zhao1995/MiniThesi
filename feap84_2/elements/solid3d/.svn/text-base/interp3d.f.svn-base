c$Id:$
      subroutine interp3d(l, xl, ndm,nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D interpolation functions for tets and bricks

c      Inputs:
c        l            - Quadrature point
c        xl(ndm,*)    - Coordinates for element nodes
c        ndm          - Space dimension of coordinates
c        nel          - Number of nodes on element

c      Outputs:
c        shp(4,nel,l) - Shape functions and derivatives
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'qudshp.h'

      integer    l, ndm,nel
      real*8     xl(ndm,*)

c     Shape functions for tetrahedron

      if(ttfl) then
        call tetshp(el3(1,l),xl,ndm,nel,jac(l),shp3(1,1,l))
        jac(l) = el3(5,l)*jac(l)

c     Shape functions for brick

      else
        call shp3d(sg3(1,l),jac(l),shp3(1,1,l),xl,ndm,nel)
        jac(l) = sg3(4,l)*jac(l)
      endif

      end
