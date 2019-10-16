c$Id:$
      subroutine strn3m(shp,xl,ul,theta,ndm,ndf,nel,nen,eps)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c     2. Compute strains at t_n                             01/06/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute mixed strain of near incompressible formulation

c      Inputs:
c        shp(4,*)      = Shape functions
c        xl(ndm,nen)   = Nodal coordinates
c        ul(ndf,nen,*) = Nodal solution parameters
c        theta(3)      = Volume change (mixed form)
c        ndm           = Spatial dimension of mesh
c        ndf           = DOF/node (max)
c        nel           = Number nodes on element (4 or 9)
c        nen           = Max nodes/element (dimension uses only)

c      Outputs:
c        eps(9,*)      = Mixed strain at point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'elcoor.h'
      include  'pconstant.h'
      include  'pmod2d.h'

      integer   ndm,      ndf,         nel,           nen,     k
      real*8    theta(3), dtheta
      real*8    shp(4,*), xl(ndm,nen), ul(ndf,nen,*), eps(9,*)

c     Compute strain tensor for constitutive equations

      do k = 1,6
        eps(k,1) = 0.0d0
        eps(k,2) = 0.0d0
      end do ! k
      do k = 1,3
        xref(k) = 0.0d0
        xcur(k) = 0.0d0
      end do ! k

      do k = 1,nel
        eps(1,1)  = eps(1,1) + shp(1,k)*ul(1,k,1)
        eps(2,1)  = eps(2,1) + shp(2,k)*ul(2,k,1)
        eps(3,1)  = eps(3,1) + shp(3,k)*ul(3,k,1)
        eps(4,1)  = eps(4,1) + shp(1,k)*ul(2,k,1)
     &                       + shp(2,k)*ul(1,k,1)
        eps(5,1)  = eps(5,1) + shp(2,k)*ul(3,k,1)
     &                       + shp(3,k)*ul(2,k,1)
        eps(6,1)  = eps(6,1) + shp(3,k)*ul(1,k,1)
     &                       + shp(1,k)*ul(3,k,1)
        eps(1,2)  = eps(1,2) + shp(1,k)*ul(1,k,2)
        eps(2,2)  = eps(2,2) + shp(2,k)*ul(2,k,2)
        eps(3,2)  = eps(3,2) + shp(3,k)*ul(3,k,2)
        eps(4,2)  = eps(4,2) + shp(1,k)*ul(2,k,2)
     &                       + shp(2,k)*ul(1,k,2)
        eps(5,2)  = eps(5,2) + shp(2,k)*ul(3,k,2)
     &                       + shp(3,k)*ul(2,k,2)
        eps(6,2)  = eps(6,2) + shp(3,k)*ul(1,k,2)
     &                       + shp(1,k)*ul(3,k,2)
        xref(1)   = xref(1)  + shp(4,k)*xl(1,k)
        xref(2)   = xref(2)  + shp(4,k)*xl(2,k)
        xref(3)   = xref(3)  + shp(4,k)*xl(3,k)
        xcur(1)   = xcur(1)  + shp(4,k)*ul(1,k,1)
        xcur(2)   = xcur(2)  + shp(4,k)*ul(2,k,1)
        xcur(3)   = xcur(3)  + shp(4,k)*ul(3,k,1)
      end do ! k

c     Set current coords

      do k = 1,3
        xcur(k) = xcur(k) + xref(k)
      end do ! k

c     Strains at t_n

      do k = 1,6
        eps(k,2) = eps(k,1) - eps(k,2)
      end do ! k

c     Correct strains and incremental strains for mixed formulation

      dtheta   = one3*(theta(1) - eps(1,1)-eps(2,1)-eps(3,1))  ! t_n+1
      eps(1,1) = eps(1,1) + dtheta
      eps(2,1) = eps(2,1) + dtheta
      eps(3,1) = eps(3,1) + dtheta

      dtheta   = theta(1) - theta(2)
      dtheta   = one3*(dtheta - eps(1,2) - eps(2,2) - eps(3,2))  ! t_n
      eps(1,2) = eps(1,2) + dtheta
      eps(2,2) = eps(2,2) + dtheta
      eps(3,2) = eps(3,2) + dtheta

      end
