c$Id:$
      subroutine sh3fener ( d, dthk, ce, cr, cx, sn, sq, sm,
     &                      ul, shp, fphi, xjw, ndf, lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute momenta and energy

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'erotas.h'
      include  'eldata.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'ptdat6.h'
      include  'tdata.h'

      integer   i         , j          , lint       , ndf
      real*8    rhoa      , rhoi       , dthk       , dot
      real*8    d    (*)  , ul  (ndf,nen,*)         , shp (4,4)
      real*8    ce   (3,4), cx  (2,4)  , cr  (3,4)  , fphi(3,4)
      real*8    sn   (3,4), sq  (2,4)  , sm  (3,4)  , aa   (3)
      real*8    pg   (3)  , vd  (3)    , w1  (3)    , xjw (4)

      save

c     Galerkin and Shear Treatment

      rhoa = d( 4)*dthk
      rhoi = d(8)*rhoa*dthk**2/12.d0

      do i = 1,lint

c       Interpolate Velocities and Accelerations

        do j = 1 , 3
          pg(j)   = shp(1,i)*fphi(j, 1)  + shp(2,i)*fphi(j, 2)
     &            + shp(3,i)*fphi(j, 3)  + shp(4,i)*fphi(j, 4)

          vd(j)   = shp(1,i)*ul  (j,1,4) + shp(2,i)*ul  (j,2,4)
     &            + shp(3,i)*ul  (j,3,4) + shp(4,i)*ul  (j,4,4)

          w1(j)   = shp(1,i)*rvel(j,1,2) + shp(2,i)*rvel(j,2,2)
     &            + shp(3,i)*rvel(j,3,2) + shp(4,i)*rvel(j,4,2)
        end do ! j

        call vecp ( pg , vd , aa )

c       Integrate Momenta

        do j = 1 , 3
          epl(j  ) = epl(j  ) + aa(j) * rhoa * xjw(i)
          epl(j+3) = epl(j+3) + w1(j) * rhoi * xjw(i)
        end do ! j

c       Integrate Energy

        epl(7) = epl(7) + 0.5d0*(dot(vd,vd,3) * rhoa
     &                         + dot(w1,w1,3) * rhoi) * xjw(i)
        epl(8) = epl(8) + 0.5d0*(dot(sn(1,i),ce(1,i),3)
     &                         + dot(sq(1,i),cx(1,i),2)
     &                         + dot(sm(1,i),cr(1,i),3)) * xjw(i)
      end do ! i

      end
