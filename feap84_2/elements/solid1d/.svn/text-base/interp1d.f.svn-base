c$Id:$
        subroutine interp1d(l, xl, ndm,nel, flag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/04/2009
c       1. Add flag to call list                            03/03/2010
c       2. Call shp1dn (natural derivatives)                14/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  1-D shape functions

c      Inputs:
c         l           - Quadrature point number
c         xl(ndm,*)   - Element nodal coordinates
c         ndm         - Space dimension of mesh
c         nel         - Number of nodes/element
c         flag        - Compute global derivatives if false

c      Outputs:
c         shp(2,*,l)  - Shape functions and first derivatives
c         jac(l)      - Jacobian of point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit     none

      include     'qudshp.h'

      logical      flag
      integer      l, ndm, nel
      real*8       xl(ndm,*)

      if(quad) then
        call shp1d(sg1(1,l),xl,shp1(1,1,l),ndm,nel,jac(l))
        jac(l) = jac(l)*sg1(2,l)
      elseif(flag) then
        call shp1dn(sg1(1,l),shp1(1,1,l),nel)
      endif

      end
