c$Id:$
      subroutine int3d(ll,lint,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Gauss quadrature for 3-d element

c      Inputs:
c         ll     - Number of points/direction

c      Outputs:
c         lint   - Total number of quadrature points
c         s(4,*) - Gauss points (1-3) and weights (4)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pconstant.h'

      integer   i,j,k,ll,lint, ig(4),jg(4)
      real*8    g, s(4,*), sw(2,5)

      save

      data      ig/-1,1,1,-1/,jg/-1,-1,1,1/

c     1 pt. quadrature

      if(ll.eq.1) then

        lint = 1
        do i = 1,3
          s(i,1) = 0.0d0
        end do ! i
        s(4,1) = 8.0d0

c     2 x 2 x 2 pt. quadrature

      elseif(ll.eq.2) then

        lint = 8
        g    = sqt13
        do i = 1,4
          s(1,i)   = ig(i)*g
          s(1,i+4) = s(1,i)
          s(2,i)   = jg(i)*g
          s(2,i+4) = s(2,i)
          s(3,i)   =  g
          s(3,i+4) = -g
          s(4,i)   = 1.d0
          s(4,i+4) = 1.d0
        end do ! i

c     Special 9 pt. quadrature

      elseif(ll.eq.-9) then

       lint = 9
       g    = sqtp6
       do i = 1,4
          s(1,i)   = ig(i)*g
          s(1,i+4) = s(1,i)
          s(2,i)   = jg(i)*g
          s(2,i+4) = s(2,i)
          s(3,i)   =  g
          s(3,i+4) = -g
          s(4,i)   = five9
          s(4,i+4) = five9
        end do ! i
        s(1,9)     =  0.d0
        s(2,9)     =  0.d0
        s(3,9)     =  0.d0
        s(4,9)     =  thty29

c     Special 4 pt. quadrature

      elseif(ll.eq.-4) then

        lint = 4
        g    = sqt13
        do i = 1,4
          s(1,i) = ig(i)*g
          s(2,i) = s(1,i)
          s(3,i) = jg(i)*g
          s(4,i) = 2.0d0
        end do ! i
        s(2,3) = -g
        s(2,4) =  g

c     ll x ll x ll pt. quadrature

      elseif(ll.le.5) then

        call int1d(ll,sw)
        lint = 0
        do k = 1,ll
          do j = 1,ll
            do i = 1,ll
              lint = lint + 1
              s(1,lint) = sw(1,i)
              s(2,lint) = sw(1,j)
              s(3,lint) = sw(1,k)
              s(4,lint) = sw(2,i)*sw(2,j)*sw(2,k)
            end do ! i
          end do ! j
        end do ! k

c     Error

      else

        write(ilg,2000) ll
        write(iow,2000) ll
        if(ior.lt.0) then
          write(*,2000) ll
        endif
        call plstop()

      endif

c     Format

2000  format('  *ERROR* INT3D: Illegal quadrature order =',i16)

      end
