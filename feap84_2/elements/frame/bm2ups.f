c$Id:$
      subroutine bm2ups(nlay,d,bz,hn,h1,nh,def,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute constitutive updates numerically

c      Inputs:
c        nlay    - Number of z-levels
c        d(*)    - Material parameters
c        bz(2,*) - Z-coordinate and B-width
c        hn(*)   - History terms at t_n
c        h1(*)   - History terms at t_n+1
c        nh      - Number of history terms per level
c        def(3)  - Axial, shear, and bending strains at t_n+alpha
c        isw     - Element control parameter

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm2com.h'

      integer   i,istrt, n,nlay,nn,nh,isw
      real*8    ta,eps,zz,d(*),bz(2,*),hn(*),h1(*),def(3)

      save

      data      ta /0.0d0/

c     Compute constitution updates using Gauss-Lobbato quadrature

      nn    = 1
      eps   = def(1) - bz(1,1)*def(3)
      istrt = nint(d(84))
      call modl1u(d,ta,eps,hn(nn),h1(nn),nh,istrt, isw)
      nn       = nn + nh
      do n = 1,nlay-1
        do i = 2,int(d(102))
          zz  = 0.5d0*((1.d0 - sl(1,i))*bz(1,n)
     &               + (1.d0 + sl(1,i))*bz(1,n+1))
          eps = def(1) - zz*def(3)
          call modl1u(d,ta,eps,hn(nn),h1(nn),nh,istrt, isw)
          nn  = nn + nh
        end do ! i
      end do ! n

      end
