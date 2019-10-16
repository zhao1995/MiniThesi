c$Id:$
      subroutine phibl(xr,xia,xib,dxi, xl, i,j,jj)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Define blended surface nodes for 3-D contact slidelines

c      Inputs:

c      Output:
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      integer    i,j,k,jj

      real*8     xr(3,181,4), xl(3,4)
      real*8     xia,xib,xi1,xi2,dxi

      xi1 = xia
      xi2 = xib
      do k = 1,3
        xl(k,1) = xr(k,i,1)*(0.5d0 - 0.5d0*xi2)
     &          + xr(k,i,3)*(0.5d0 + 0.5d0*xi2)
     &          + xr(k,j,4)*(0.5d0 - 0.5d0*xi1)
     &          + xr(k,j,2)*(0.5d0 + 0.5d0*xi1)
     &          - xr(k, 1,1)*(0.5d0 - 0.5d0*xi1)*(0.5d0 - 0.5d0*xi2)
     &          - xr(k,jj,1)*(0.5d0 + 0.5d0*xi1)*(0.5d0 - 0.5d0*xi2)
     &          - xr(k, 1,3)*(0.5d0 - 0.5d0*xi1)*(0.5d0 + 0.5d0*xi2)
     &          - xr(k,jj,3)*(0.5d0 + 0.5d0*xi1)*(0.5d0 + 0.5d0*xi2)
      end do

      xi1 = xia + dxi
      xi2 = xib
      do k = 1,3
        xl(k,2) = xr(k,i+1,1)*(0.5d0 - 0.5d0*xi2)
     &          + xr(k,i+1,3)*(0.5d0 + 0.5d0*xi2)
     &          + xr(k,j,4)*(0.5d0 - 0.5d0*xi1)
     &          + xr(k,j,2)*(0.5d0 + 0.5d0*xi1)
     &          - xr(k, 1,1)*(0.5d0 - 0.5d0*xi1)*(0.5d0 - 0.5d0*xi2)
     &          - xr(k,jj,1)*(0.5d0 + 0.5d0*xi1)*(0.5d0 - 0.5d0*xi2)
     &          - xr(k, 1,3)*(0.5d0 - 0.5d0*xi1)*(0.5d0 + 0.5d0*xi2)
     &          - xr(k,jj,3)*(0.5d0 + 0.5d0*xi1)*(0.5d0 + 0.5d0*xi2)
      end do

      xi1 = xia + dxi
      xi2 = xib + dxi
      do k = 1,3
        xl(k,3) = xr(k,i+1,1)*(0.5d0 - 0.5d0*xi2)
     &          + xr(k,i+1,3)*(0.5d0 + 0.5d0*xi2)
     &          + xr(k,j+1,4)*(0.5d0 - 0.5d0*xi1)
     &          + xr(k,j+1,2)*(0.5d0 + 0.5d0*xi1)
     &          - xr(k, 1,1)*(0.5d0 - 0.5d0*xi1)*(0.5d0 - 0.5d0*xi2)
     &          - xr(k,jj,1)*(0.5d0 + 0.5d0*xi1)*(0.5d0 - 0.5d0*xi2)
     &          - xr(k, 1,3)*(0.5d0 - 0.5d0*xi1)*(0.5d0 + 0.5d0*xi2)
     &          - xr(k,jj,3)*(0.5d0 + 0.5d0*xi1)*(0.5d0 + 0.5d0*xi2)
      end do

      xi1 = xia
      xi2 = xib + dxi
      do k = 1,3
        xl(k,4) = xr(k,i,1)*(0.5d0 - 0.5d0*xi2)
     &          + xr(k,i,3)*(0.5d0 + 0.5d0*xi2)
     &          + xr(k,j+1,4)*(0.5d0 - 0.5d0*xi1)
     &          + xr(k,j+1,2)*(0.5d0 + 0.5d0*xi1)
     &          - xr(k, 1,1)*(0.5d0 - 0.5d0*xi1)*(0.5d0 - 0.5d0*xi2)
     &          - xr(k,jj,1)*(0.5d0 + 0.5d0*xi1)*(0.5d0 - 0.5d0*xi2)
     &          - xr(k, 1,3)*(0.5d0 - 0.5d0*xi1)*(0.5d0 + 0.5d0*xi2)
     &          - xr(k,jj,3)*(0.5d0 + 0.5d0*xi1)*(0.5d0 + 0.5d0*xi2)
      end do

      end
