c$Id:$
      subroutine pshp9(xi,shp,d1shp,d2shp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute 9-node shape functions.

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    xi(2), shp(9),d1shp(2,9),d2shp(3,9)
      real*8    sh1m,sh2m,sh1p,sh2p,sh1c,sh2c
      real*8    dn1m,dn2m,dn1p,dn2p,dn1c,dn2c, ddn1,ddn2,ddnc

      save

c     Constant parameters for 9-node shape functions

      sh1m = 0.5d0*(xi(1)*xi(1) - xi(1))
      sh2m = 0.5d0*(xi(2)*xi(2) - xi(2))
      sh1p = 0.5d0*(xi(1)*xi(1) + xi(1))
      sh2p = 0.5d0*(xi(2)*xi(2) + xi(2))
      sh1c = 1.0d0 - xi(1)*xi(1)
      sh2c = 1.0d0 - xi(2)*xi(2)

      dn1m = xi(1) - 0.5d0
      dn2m = xi(2) - 0.5d0
      dn1p = xi(1) + 0.5d0
      dn2p = xi(2) + 0.5d0
      dn1c = -2.0d0*xi(1)
      dn2c = -2.0d0*xi(2)

      ddn1 =  1.0d0
      ddn2 =  1.0d0
      ddnc = -2.0d0

c     Form shape functions

      shp(1)     = sh1m*sh2m
      shp(2)     = sh1p*sh2m
      shp(3)     = sh1p*sh2p
      shp(4)     = sh1m*sh2p
      shp(5)     = sh1c*sh2m
      shp(6)     = sh1p*sh2c
      shp(7)     = sh1c*sh2p
      shp(8)     = sh1m*sh2c
      shp(9)     = sh1c*sh2c

c     Form first derivatives of shape functions

      d1shp(1,1) = dn1m*sh2m
      d1shp(1,2) = dn1p*sh2m
      d1shp(1,3) = dn1p*sh2p
      d1shp(1,4) = dn1m*sh2p
      d1shp(1,5) = dn1c*sh2m
      d1shp(1,6) = dn1p*sh2c
      d1shp(1,7) = dn1c*sh2p
      d1shp(1,8) = dn1m*sh2c
      d1shp(1,9) = dn1c*sh2c

      d1shp(2,1) = sh1m*dn2m
      d1shp(2,2) = sh1p*dn2m
      d1shp(2,3) = sh1p*dn2p
      d1shp(2,4) = sh1m*dn2p
      d1shp(2,5) = sh1c*dn2m
      d1shp(2,6) = sh1p*dn2c
      d1shp(2,7) = sh1c*dn2p
      d1shp(2,8) = sh1m*dn2c
      d1shp(2,9) = sh1c*dn2c

c     Form second derivatives of shape functions

      d2shp(1,1) = ddn1*sh2m
      d2shp(1,2) = ddn1*sh2m
      d2shp(1,3) = ddn1*sh2p
      d2shp(1,4) = ddn1*sh2p
      d2shp(1,5) = ddnc*sh2m
      d2shp(1,6) = ddn1*sh2c
      d2shp(1,7) = ddnc*sh2p
      d2shp(1,8) = ddn1*sh2c
      d2shp(1,9) = ddnc*sh2c

      d2shp(2,1) = dn1m*dn2m
      d2shp(2,2) = dn1p*dn2m
      d2shp(2,3) = dn1p*dn2p
      d2shp(2,4) = dn1m*dn2p
      d2shp(2,5) = dn1c*dn2m
      d2shp(2,6) = dn1p*dn2c
      d2shp(2,7) = dn1c*dn2p
      d2shp(2,8) = dn1m*dn2c
      d2shp(2,9) = dn1c*dn2c

      d2shp(3,1) = sh1m*ddn2
      d2shp(3,2) = sh1p*ddn2
      d2shp(3,3) = sh1p*ddn2
      d2shp(3,4) = sh1m*ddn2
      d2shp(3,5) = sh1c*ddn2
      d2shp(3,6) = sh1p*ddnc
      d2shp(3,7) = sh1c*ddn2
      d2shp(3,8) = sh1m*ddnc
      d2shp(3,9) = sh1c*ddnc

      end
