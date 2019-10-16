c$Id:$
      subroutine pgauss(l,lint,r,z,w)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form Gauss points and weights for two dimensions

c      Inputs:
c         l       - Number of points/direction

c      Outputs:
c         lint    - Total number of points
c         r(*)    - 1-direction Gauss point
c         z(*)    - 2-direction Gauss point
c         w(*)    - Gauss weight
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'

      integer   i,j,k,l,lint, lr(9),lz(9),lw(9)
      real*8    g,h, r(*),z(*),w(*),sg(5),wg(5)

      save

      data      lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
      data      lw/4*25,4*40,64/

c     Set total number of points

      lint = l*l

c     5 pt. integration

      if(l.eq.0) then

        lint = 5
        g    = sqrt(0.6d0)
        do i = 1,4
          r(i) = g*lr(i)
          z(i) = g*lz(i)
          w(i) = 5.d0/9.d0
        end do ! i

        r(5) = 0.0d0
        z(5) = 0.0d0
        w(5) = 16.d0/9.d0

c     1x1 integration

      elseif(l.eq.1) then
        r(1) = 0.d0
        z(1) = 0.d0
        if(nel.eq.3) z(1) = -1./3.0d0
        w(1) = 4.d0

c     2x2 integration

      elseif(l.eq.2) then
        g = 1.d0/sqrt(3.d0)
        do i = 1,4
          r(i) = g*lr(i)
          z(i) = g*lz(i)
          w(i) = 1.d0
        end do ! i

c     3x3 integration

      elseif(l.eq.3) then
        g = sqrt(0.6d0)
        h = 1.d0/81.d0
        do i = 1,9
          r(i) = g*lr(i)
          z(i) = g*lz(i)
          w(i) = h*lw(i)
        end do ! i

c     4x4 integration

      elseif(l.eq.4) then
        g = sqrt(4.8d0)
        h = sqrt(30.0d0)/36.0d0
        sg(1) = sqrt((3.d0+g)/7.d0)
        sg(4) = - sg(1)
        sg(2) = sqrt((3.d0-g)/7.d0)
        sg(3) = -sg(2)
        wg(1) = 0.5d0 - h
        wg(2) = 0.5d0 + h
        wg(3) = 0.5d0 + h
        wg(4) = 0.5d0 - h
        i = 0
        do j = 1,4
          do k = 1,4
            i = i + 1
            r(i) = sg(k)
            z(i) = sg(j)
            w(i) = wg(j)*wg(k)
          end do ! k
        end do ! i

c     5x5 integration

      elseif(l.eq.5) then

        g     = sqrt(1120.d0)
        sg(1) =  sqrt((70.d0 + g)/126.d0)
        sg(2) =  sqrt((70.d0 - g)/126.d0)
        sg(3) =  0.0d0
        sg(4) = -sg(2)
        sg(5) = -sg(1)

        wg(1) =  (21.d0*g + 117.6d0)/(g*(70.d0 + g))
        wg(2) =  (21.d0*g - 117.6d0)/(g*(70.d0 - g))
        wg(3) =  2.d0*(1.d0 - wg(1) - wg(2))
        wg(4) =  wg(2)
        wg(5) =  wg(1)

        i = 0
        do j = 1,5
          do k = 1,5
            i = i + 1
            r(i) = sg(k)
            z(i) = sg(j)
            w(i) = wg(j)*wg(k)
          end do ! k
        end do ! j

      endif

      end
