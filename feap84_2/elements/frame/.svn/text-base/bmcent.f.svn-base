c$Id:$
      subroutine bmcent(n,d)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute geometric properties for centroid

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   i,n
      real*8    a,az,bi,bp,zi,zp,b1,b2,z1,zbar, d(2,*)

      save

c     Integrate for area and first moment of area

      a  = 0.0d0
      az = 0.0d0
      do i = 1,n-1
        zi = d(1,i)
        bi = d(2,i)
        zp = d(1,i+1)
        bp = d(2,i+1)
        b1 = bi*zp - bp*zi
        b2 = bp - bi
        z1 = 0.5d0*(zp + zi)
        a  = a  + b1    + b2*z1
        az = az + b1*z1 + b2*(zp*(zp + zi) + zi*zi)*one3
      end do ! i

c     Centroid location

      zbar = az/a

c     Shift z-coordinates by centroid value

      do i = 1,n
        d(1,i) = d(1,i) - zbar
      end do ! i

      end
