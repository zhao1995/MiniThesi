c$Id:$
      subroutine sengy(ul,shp, dvol,dmas,lfac,cfac, is,ns, mom)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change shp(4, to shp(ns,                         19/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'cdata.h'
      include   'eldata.h'
      include   'elengy.h'
      include   'ptdat6.h'
      include   'sdata.h'

      integer    is,ns, i,j
      real*8     vnode,vl,vsqr, dvol,dmas, lfac,cfac, cons
      real*8     ul(ndf,nen,*),shp(ns,*),mom(3,*)

c     Compute sum of velocity squared at point

      vsqr = 0.0d0
      do i = 1,is
        vl = 0.0d0
        do j = 1,nel
          vl = vl + ul(i,j,4)*shp(ns,j)
        end do ! j
        vsqr = vsqr + vl*vl
        vl   = vl*cfac*dmas                  ! Set consistent mass part
        do j = 1,nel
          mom(i,j) = mom(i,j) + shp(ns,j)*vl
        end do ! j
      end do ! i

c     Lumped mass contribution

      vnode = 0.0d0
      do j = 1,nel
        vl = 0.0d0
        do i = 1,is
          vl = vl + ul(i,j,4)**2
        end do ! i
        vnode = vnode + vl*shp(ns,j)
        cons  = shp(ns,j)*lfac*dmas
        do i = 1,is
          mom(i,j) = mom(i,j) + ul(i,j,4)*cons
        end do ! i
      end do ! j

c     Accumulate kinetic energy

      epl(7) = epl(7) + 0.5d0*(lfac*vnode + cfac*vsqr)*dmas

c     Accumulate stored energy

      epl(8) = epl(8) + estore*dvol

      end
