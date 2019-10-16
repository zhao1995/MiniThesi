c$Id:$
      subroutine strn1m(shp,xl,ul,theta,irad,ndm,ndf,nel,nen,eps)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute mixed strain of near incompressible formulation

c      Inputs:
c        shp(2,*)      = Shape functions
c        xl(ndm,nen)   = Nodal coordinates
c        ul(ndf,nen,*) = Nodal solution parameters
c        theta         = Volume change (mixed form)
c        irad          = Inverse radius (or zero for plane).
c        ndm           = Spatial dimension of mesh
c        ndf           = DOF/node (max)
c        nel           = Number nodes on element
c        nen           = Max nodes/element (dimension uses only)

c      Outputs:
c        eps(*)        = Mixed strain at point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'elcoor.h'
      include  'pconstant.h'

      integer   ndm,      ndf,         nel,         nen,     k
      real*8    irad,     theta,       dtheta
      real*8    shp(2,*), xl(ndm,nen), ul(ndf,nen), eps(*)

c     Compute strain tensor for constitutive equations

      do k = 1,3
        eps(k)  = 0.0d0
        xref(k) = 0.0d0
        xcur(k) = 0.0d0
      end do ! k

      do k = 1,nel
        eps(1)  = eps(1)  + shp(1,k)*ul(1,k)
        eps(3)  = eps(3)  + shp(2,k)*ul(1,k)
        xref(1) = xref(1) + shp(2,k)*xl(1,k)
        xcur(1) = xcur(1) + shp(2,k)*ul(1,k)
      end do ! k

c     Set current coords

      xcur(1) = xcur(1) + xref(1)

c     Modification for plane/axisymmetry

      eps(3) = eps(3)*irad

c     Correct strains and incremental strains for mixed formulation

      dtheta = one3*(theta - eps(1) - eps(2) - eps(3))
      eps(1) = eps(1) + dtheta
      eps(2) = eps(2) + dtheta
      eps(3) = eps(3) + dtheta

      end
