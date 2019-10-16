c$Id:$
      subroutine int2d(l,lint,sg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'one3' from 'pconstant.h' instead of 'third'   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form Gauss points and weights for two dimensions

c      Inputs:
c         l       - Number of points/direction

c      Outputs:
c         lint    - Total number of points
c         sg(3,*) - Array of points and weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'iofile.h'
      include  'pconstant.h'

      integer   i,j,k,l,lint, lr(9),lz(9),lw(9)
      real*8    g,h, sg(3,*),ss(5),ww(5)

      save

      data      lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
      data      lw/4*25,4*40,64/

c     Set number of total points

      lint = l*l

c     5 pt. integration

      if(l.eq.0) then

        lint = 5
        g    = sqtp6                          ! sqrt(0.6)
        do i = 1,4
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = five9
        end do ! i

        sg(1,5) = 0.0d0
        sg(2,5) = 0.0d0
        sg(3,5) = 2.8d0*eight9                ! 16/9

c     1x1 integration

      elseif(l.eq.1) then
        sg(1,1) = 0.d0
        sg(2,1) = 0.d0
        if(nel.eq.3) sg(2,1) = -one3
        sg(3,1) = 4.d0

c     2x2 integration

      elseif(l.eq.2) then
        g = sqt13                              ! sqrt(1/3)
        do i = 1,4
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = 1.d0
        end do ! i

c     3x3 integration

      elseif(l.eq.3) then
        g = sqtp6                              ! sqrt(0.6)
        h = 1.d0
        h = h/81.d0
        do i = 1,9
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = h*lw(i)
        end do ! i

c     4x4 integration

      elseif(l.eq.4) then
        g     = sqt48                         ! sqrt(4.8)
        h     = one3/g
        ss(1) = sqrt((3.d0+g)/7.d0)
        ss(4) = -ss(1)
        ss(2) = sqrt((3.d0-g)/7.d0)
        ss(3) = -ss(2)
        ww(1) = 0.5d0 - h
        ww(2) = 0.5d0 + h
        ww(3) = 0.5d0 + h
        ww(4) = 0.5d0 - h
        i = 0
        do j = 1,4
          do k = 1,4
            i = i + 1
            sg(1,i) = ss(k)
            sg(2,i) = ss(j)
            sg(3,i) = ww(j)*ww(k)
          end do ! k
        end do ! i

c     5x5 integration

      elseif(l.eq.5) then

        g     =  1120.d0
        g     =  sqrt(g)
        ss(1) =  sqrt((70.d0 + g)/126.d0)
        ss(2) =  sqrt((70.d0 - g)/126.d0)
        ss(3) =  0.0d0
        ss(4) = -ss(2)
        ss(5) = -ss(1)

        ww(1) =  (21.d0*g + 117.6d0)/(g*(70.d0 + g))
        ww(2) =  (21.d0*g - 117.6d0)/(g*(70.d0 - g))
        ww(3) =  2.d0*(1.d0 - ww(1) - ww(2))
        ww(4) =  ww(2)
        ww(5) =  ww(1)

        i = 0
        do j = 1,5
          do k = 1,5
            i = i + 1
            sg(1,i) = ss(k)
            sg(2,i) = ss(j)
            sg(3,i) = ww(j)*ww(k)
          end do ! k
        end do ! j

c     Error

      else

        write(ilg,2000) l
        write(iow,2000) l
        if(ior.lt.0) then
          write(*,2000) l
        endif
        call plstop()

      endif

c     Format

2000  format(' *ERROR* INT2D: Illegal quadrature order =',i16)

      end
