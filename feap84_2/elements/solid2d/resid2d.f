c$Id:$
      subroutine resid2d(cfac,lfac,xsj,xsj0,xx,yy,shp,sig,d,vl,al,p,ndf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add array 'bdy' to control body loading          20/07/2007
c          Add radial body loading option
c       2. Add computation of 'bdy' for other cases         17/04/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plane and axisymmetric residual routine

c      Inputs:

c      Outputs:
c         p(ndf,*)  - Element residual
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'eltran.h'
      include  'pmod2d.h'

      integer   ndf, j
      real*8    aj1,aj2,aj3,aj0,rr,lfac,cfac,xsj,xsj0,xx,yy
      real*8    d(*),vl(ndf,*),al(ndf,*),p(ndf,*)
      real*8    shp(3,*),sig(*),vc(3),ac(3),body(3),bdy(3)

      save

c     Compute body force values

      call sbodyf(d, body)
      if(stype.lt.3 .and. nint(d(69)).eq.5) then
        rr = sqrt(xx*xx + yy*yy)
        if(rr.gt.0.0d0) then
          bdy(1) = (xx*body(1) - yy*body(2))/rr
          bdy(2) = (yy*body(1) + xx*body(2))/rr
        else
          do j = 1,3
            bdy(j) = body(j)
          end do ! j
        endif
      else
        do j = 1,3
          bdy(j) = body(j)
        end do ! j
      endif

c     Compute accelerations

      ac(1) = 0.0d0
      ac(2) = 0.0d0
      do j = 1,nel
        ac(1) = ac(1) + shp(3,j)*al(1,j)
        ac(2) = ac(2) + shp(3,j)*al(2,j)
      end do ! j

      rr    = d(4)
      ac(1) = rr*ac(1)*cfac
      ac(2) = rr*ac(2)*cfac

c     For Rayleigh Mass Damping: Compute velocity

      if(d(77).ne.0.0d0) then
        vc(1) = 0.0d0
        vc(2) = 0.0d0
        do j = 1,nel
          vc(1) = vc(1) + shp(3,j)*vl(1,j)
          vc(2) = vc(2) + shp(3,j)*vl(2,j)
        end do ! j
        vc(1) = rr*vc(1)*cfac*d(77)
        vc(2) = rr*vc(2)*cfac*d(77)

        aj0 = lfac*d(77)*rr
        do j = 1,nel
          p(1,j) = p(1,j) - (vc(1) + aj0*vl(1,j))*shp(3,j)*xsj
          p(2,j) = p(2,j) - (vc(2) + aj0*vl(2,j))*shp(3,j)*xsj
        end do ! j

      endif

c     Loop over rows

      do j = 1,nel
        aj1 = shp(1,j)*xsj
        aj2 = shp(2,j)*xsj
        aj3 = shp(3,j)*xsj0
        aj0 = lfac*rr

c       Compute gravity, thermal, inertia, and stress contributions

        p(1,j) = p(1,j) + (bdy(1) - ac(1) - aj0*al(1,j))*shp(3,j)*xsj
     &                  -  aj1*sig(1) - aj2*sig(4) - aj3*sig(3)
        p(2,j) = p(2,j) + (bdy(2) - ac(2) - aj0*al(2,j))*shp(3,j)*xsj
     &                  -  aj1*sig(4) - aj2*sig(2)

      end do ! j

c     Axisymmetry with torsion

      if(stype.eq.8) then
        ac(3) = 0.0d0
        vc(3) = 0.0d0
        do j = 1,nel
          ac(3) = ac(3) + shp(3,j)*al(3,j)
          vc(3) = vc(3) + shp(3,j)*vl(3,j)
        end do ! j
        ac(3)   = rr*ac(3)*cfac
        vc(3)   = rr*vc(3)*cfac
        aj0     = rr*lfac
        do j = 1,nel
          p(3,j) = p(3,j) + ((bdy(3) - (ac(3) + aj0*al(3,j))
     &                    -  d(77)*(vc(3) + aj0*vl(3,j)))*shp(3,j)
     &                    -  shp(2,j)*sig(5) -  shp(1,j)*sig(6))*xsj
     &                    +  shp(3,j)*sig(6)*xsj0
        end do ! j
      endif

c     Rotational forces: ro * r * omega^2

      if(xsj0.eq.0.0d0) then
        do j = 1,nel
          aj0    = shp(3,j)*d(4)*d(65)**2*xsj
          p(1,j) = p(1,j) + aj0*xx
          p(2,j) = p(2,j) + aj0*yy
        end do ! j
      else
        do j = 1,nel
          p(1,j) = p(1,j) + shp(3,j)*d(4)*xx*d(65)**2*xsj
        end do ! j
      endif

c     Axisymmetric patch loading

      if(stype.eq.3) then
        if(max(abs(d(197)),abs(d(198))).gt.0.0d0) then
          do j = 1,nel
            p(1,j) = p(1,j) + shp(3,j)*d(197)/(xx**2)*xsj
            p(2,j) = p(2,j) + shp(3,j)*d(198)*xsj
          end do ! i
        endif
      endif

      end
