c$Id:$
      subroutine strn2m(shp,xl,ul,theta,irad,ndm,ndf,nel,nen,eps)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Compute strains at t_n                           01/06/2009
c          Add common 'pconstant.h'
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute mixed strain of near incompressible formulation

c      Inputs:
c        shp(3,nen,*)  = Shape functions
c        xl(ndm,nen)   = Nodal coordinates
c        ul(ndf,nen,*) = Nodal solution parameters
c        theta(2)      = Volume change (mixed form)
c        irad          = Inverse radius (or zero for plane).
c        ndm           = Spatial dimension of mesh
c        ndf           = DOF/node (max)
c        nel           = Number nodes on element (4 or 9)
c        nen           = Max nodes/element (dimension uses only)

c      Outputs:
c        eps(9,3)      = Mixed strain at point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'elcoor.h'
      include  'pmod2d.h'
      include  'pconstant.h'

      integer   ndm,      ndf,         nel,         nen,     k
      real*8    irad,     theta(2),    dtheta
      real*8    shp(3,*), xl(ndm,nen), ul(ndf,nen,*), eps(9,*)

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
        eps(1,1)  = eps(1,1)  + shp(1,k)*ul(1,k,1)
        eps(2,1)  = eps(2,1)  + shp(2,k)*ul(2,k,1)
        eps(3,1)  = eps(3,1)  + shp(3,k)*ul(1,k,1)
        eps(4,1)  = eps(4,1)  + shp(1,k)*ul(2,k,1)
     &                        + shp(2,k)*ul(1,k,1)
        eps(1,2)  = eps(1,2)  + shp(1,k)*ul(1,k,2)
        eps(2,2)  = eps(2,2)  + shp(2,k)*ul(2,k,2)
        eps(3,2)  = eps(3,2)  + shp(3,k)*ul(1,k,2)
        eps(4,2)  = eps(4,2)  + shp(1,k)*ul(2,k,2)
     &                        + shp(2,k)*ul(1,k,2)
        xref(1)   = xref(1)   + shp(3,k)*xl(1,k)
        xref(2)   = xref(2)   + shp(3,k)*xl(2,k)
        xcur(1)   = xcur(1)   + shp(3,k)*ul(1,k,1)
        xcur(2)   = xcur(2)   + shp(3,k)*ul(2,k,1)
      end do ! k

c     Set current coords

      xcur(1) = xcur(1) + xref(1)
      xcur(2) = xcur(2) + xref(2)

c     Modification for plane/axisymmetry

      eps(3,1) = eps(3,1)*irad
      eps(3,2) = eps(3,2)*irad

c     Compute strains at t_n

      do k = 1,4
        eps(k,2) = eps(k,1) - eps(k,2)
      end do ! k

c     Correct strains and incremental strains for mixed formulation

      dtheta   = one3*(theta(1) - eps(1,1) - eps(2,1) - eps(3,1))
      eps(1,1) = eps(1,1) + dtheta
      eps(2,1) = eps(2,1) + dtheta
      eps(3,1) = eps(3,1) + dtheta

      dtheta   = theta(1) - theta(2)
      dtheta   = one3*(dtheta - eps(1,2) - eps(2,2) - eps(3,2))
      eps(1,2) = eps(1,2) + dtheta
      eps(2,2) = eps(2,2) + dtheta
      eps(3,2) = eps(3,2) + dtheta

c     Torsion case

      if(stype.eq.8) then
        do k = 1,nel
          eps(5,1) = eps(5,1) + shp(2,k)*ul(3,k,1)
          eps(6,1) = eps(6,1) + (shp(1,k) - irad*shp(3,k))*ul(3,k,1)
          eps(5,2) = eps(5,2) + shp(2,k)*ul(3,k,2)
          eps(6,2) = eps(6,2) + (shp(1,k) - irad*shp(3,k))*ul(3,k,2)
        end do ! K
        eps(5,2) = eps(5,1) - eps(5,2)
        eps(6,2) = eps(6,1) - eps(6,2)
      endif

      end
