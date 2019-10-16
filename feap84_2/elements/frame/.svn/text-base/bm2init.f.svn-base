c$Id:$
      subroutine bm2init (d,hn,h1,nh)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Initial property inputs

c     Inputs:
c        d(*)      - Material parameters
c        hn(*)     - History terms at t_n
c        h1(*)     - History terms at t_n+1
c        nh        - Number of history terms per level

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nh,nlay, istrt, n,nn,i,ii
      real*8    d(*),hn(*),h1(*)
      real*8    ta,eps(2),sig(2),dd(2)

      save

      data     ta /0.0d0/, eps /0.0d0, 0.0d0/

c     Thickness quadrature model

      nlay = int(d(101))
      if(nlay.gt.0) then

c       Initialize constitution for Gauss-Lobbato quadrature

        nn     = 1
        ii     = 1
        istrt  = nint(d(84))
        call modl1d(d,ta,eps,hn(nn),h1(nn),nh,ii,istrt, dd,sig, 14)

        do n = 1,nlay-1

c         Add remaining points

          do i = 2,nint(d(102))

c           Increment counters

            nn  = nn + nh
            ii  = ii + 1

c           Initialize constitution

            call modl1d(d,ta,eps,hn(nn),h1(nn),nh,ii,istrt, dd,sig, 14)

          end do ! i
        end do ! n

      endif

      end
