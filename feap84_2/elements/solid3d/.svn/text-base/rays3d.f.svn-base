c$Id:$
      subroutine rays3d(d,shp,shpbar,sig,dr,vl,ndf,nel,mixed)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Stiffness proportional Rayleigh damping residual

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      logical   mixed
      integer   ndf,nel, i,j
      real*8    d(*),shp(4,*),shpbar(3,*),sig(*),dr(6,6),vl(ndf,*)
      real*8    eps(6), theta,dtheta

c     Compute strain rate terms

      do j = 1,6
        eps(j) = 0.0d0
      end do ! j
      do j = 1,nel
        eps(1) = eps(1) + shp(1,j)*vl(1,j)
        eps(2) = eps(2) + shp(2,j)*vl(2,j)
        eps(3) = eps(3) + shp(3,j)*vl(3,j)
        eps(4) = eps(4) + shp(2,j)*vl(1,j) + shp(1,j)*vl(2,j)
        eps(5) = eps(5) + shp(3,j)*vl(2,j) + shp(2,j)*vl(3,j)
        eps(6) = eps(6) + shp(1,j)*vl(3,j) + shp(3,j)*vl(1,j)
      end do ! j

c     Modify if mixed

      if(mixed) then

c       Mixed volume change

        theta  = 0.0d0
        do j = 1,nel
          theta  = theta  + shpbar(1,j)*vl(1,j)
     &                    + shpbar(2,j)*vl(2,j)
     &                    + shpbar(3,j)*vl(3,j)
        end do ! j

c       Mixed strains

        dtheta = one3*(theta - eps(1) - eps(2) - eps(3))
        eps(1) = eps(1) + dtheta
        eps(2) = eps(2) + dtheta
        eps(3) = eps(3) + dtheta

      endif

c     Compute stress modification due to Rayleigh damping

      do j = 1,6
        eps(j) = eps(j)*d(78)
        do i = 1,6
          sig(i) = sig(i) + dr(i,j)*eps(j)
        end do ! i
      end do ! j

      end
