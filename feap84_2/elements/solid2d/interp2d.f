c$Id:$
      subroutine interp2d(l, xl,ix, ndm,nel, flag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    11/11/2008
c       1. Add B-spline option for 6-node triangles         01/05/2010
c       2. Add 'shptri' for natural derivative only option  11/12/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Interpolation functions for 2-D elements

c      Inputs:
c         l            - Quadrature point
c         xl(ndm,*)    - Nodal coordinates
c         ix(*)        - Global nodal connections
c         ndm          - Mesh coordinate dimension
c         nel          - Number of element nodes
c         flag         - Global derivatives if .false.

c      Outputs: Through common block /qudshp*/
c         shp(3,16,l)  - Shape functions
c         jac(l)       - Jacobian
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'qudshp.h'

      logical    flag
      integer    l, ndm,nel, ix(*)
      real*8     xl(ndm,*)

      if(quad) then         ! Quadrilateral element
        call shp2d(sg2(1,l),xl,shp2(1,1,l),jac(l),ndm,nel,ix,flag)
        jac(l) = jac(l)*sg2(3,l)
      elseif(bsplfl) then   ! 6-node B-spline triangle
        call trispl(el2(1,l),xl,ndm, jac(l),shp2(1,1,l))
        jac(l) = jac(l)*el2(4,l)
      else                  ! Triangular element
        if(flag) then
          call shptri(el2(1,l), nel, shp2(1,1,l))
          sg2(3,l) = 0.5d0*el2(4,l) ! For proper area of cross products
        else
          call trishp(el2(1,l),xl,ndm,nel-4,jac(l),shp2(1,1,l))
          jac(l) = jac(l)*el2(4,l)
        endif
      endif

      end
