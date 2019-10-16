c$Id:$
      subroutine bm2rho(nlay,d,bz, rhoa,rhoi)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute resultant inertia properties

c     Inputs:
c        nlay    - Number of z-levels
c        d(*)    - Material parameters
c        bz(2,*) - Z-coordinate and B-width

c     Outputs:
c        rhoa    - integral * rho * b(z) * dz
c        rhoi    - integral * rho * b(z) * z*z * dz
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm2com.h'

      integer   i,n,nlay
      real*8    rhoa,rhoi, zz,df,dz, d(*),bz(2,*),brho(50)

      save

c     Compute density times width at each level

      do n = 1,nlay
        brho(n) = bz(2,n)*d(4)
      end do ! n

c     Initialize values

      rhoa = 0.0d0
      rhoi = 0.0d0

c     Integrate over each layer

      do n = 1,nlay-1
        dz = 0.5d0*(bz(1,n+1) - bz(1,n))
        if(dz .gt. 0.0d0) then
          do i = 1,int(d(102))
            zz   = 0.5d0*((1.d0 - sl(1,i))*bz(1,n)
     &                  + (1.d0 + sl(1,i))*bz(1,n+1))
            df   = 0.5d0*((1.d0 - sl(1,i))*brho(n)
     &                  + (1.d0 + sl(1,i))*brho(n+1))*dz*sl(2,i)

            rhoa =  rhoa + df
            rhoi =  rhoi + df*zz*zz
          end do ! i
        endif
      end do ! n

      end
