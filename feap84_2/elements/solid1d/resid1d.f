c$Id:$
      subroutine resid1d(cfac,lfac,xsj,xsj0,shp,sig,d,vl,al,p,ndf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plane and axisymmetric residual routine

c      Inputs:

c      Outputs:
c         p(ndf,*)  - Element residual
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'eltran.h'

      integer   ndf, j
      real*8    aj1,aj2,aj0,rr,lfac,cfac,xsj,xsj0, vc,ac
      real*8    d(*),vl(ndf,*),al(ndf,*),p(ndf,*),shp(2,*),sig(*),bf(3)

      save

c     Compute body force values

      call sbodyf(d, bf)

c     Compute accelerations

      ac = 0.0d0
      do j = 1,nel
        ac = ac + shp(2,j)*al(1,j)
      end do ! j
      rr = d(4)
      ac = rr*ac*cfac

c     For Rayleigh Mass Damping: Compute velocity

      if(d(77).ne.0.0d0) then
        vc = 0.0d0
        do j = 1,nel
          vc = vc + shp(2,j)*vl(1,j)
        end do ! j
        vc = rr*vc*cfac*d(77)

        do j = 1,nel
          aj0 = lfac*d(77)*rr
          p(1,j) = p(1,j) - (vc + aj0*vl(1,j))*shp(2,j)*xsj
        end do ! j

      endif

c     Loop over rows

      do j = 1,nel
        aj1 = shp(1,j)*xsj
        aj2 = shp(2,j)*xsj0
        aj0 = lfac*rr

c       Compute gravity, thermal, inertia, and stress contributions

        p(1,j) = p(1,j) + (bf(1) - ac - aj0*al(1,j))*shp(2,j)*xsj
     &                  - aj1*sig(1) - aj2*sig(3)
      end do ! j

      end
