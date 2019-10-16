c$Id:$
      subroutine bmat1d(c,r,shp,g,bbar)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sets B-bar matrix for 1-d problems

c      Inputs:
c         c         - Constant for plane = 0; for axisymm = 1
c         r         - Radius for axisymmetrix (= 1 for plane)
c         shp(2)    - Shape function and derivatives
c         g         - b-bar integrals

c      Outputs:
c         bbar(4)   - B-bar matrix for a node
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      real*8    c,r,bb1,sh2, shp(2),g,bbar(4)

c     Mixed modification to form B-bar

      sh2 = c*shp(2)/r
      bb1 = (g - shp(1) - sh2)*one3

c     B-bar matrix for plane and axisymmetric problems

      bbar(1) = bb1 + shp(1)
      bbar(2) = bb1
      bbar(3) = bb1 + sh2

      end
