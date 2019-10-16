c$Id:$
      subroutine tint2dn(l,lint,el)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Replace constant by 'one3'                         14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set quadrature points & weights for triangular element

c      Inputs:
c         l       - Number of nodes on element

c      Outputs:
c         lint    - Total number of points
c         el(4,*) - Area coordinate points and weights for quadrature
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pconstant.h'

      integer   i,j,l, lint
      real*8    el(4,*), gg

      save

c     3-point nodal integration

      if(l.eq.3) then
        do j = 1,3
          do i = 1,3
            el(i,j) = 0.0d0
          end do ! i
          el(j,j) = 1.d0
          el(4,j) = one3
        end do ! j

        lint    = 3

c     6-point nodal integration

      elseif(l.eq.6) then
        gg = 1.d0/15.d0
        do j = 1,3
          do i = 1,3
            el(i,j) = 0.0d0
          end do ! i
          el(j,j) = 1.d0
          el(4,j) = gg
        end do ! j
        el(1,4) = 0.5d0
        el(2,4) = 0.5d0
        el(3,4) = 0.0d0
        el(4,4) = 4.d0*gg

        el(1,5) = 0.0d0
        el(2,5) = 0.5d0
        el(3,5) = 0.5d0
        el(4,5) = el(4,4)

        el(1,6) = 0.5d0
        el(2,6) = 0.0d0
        el(3,6) = 0.5d0
        el(4,6) = el(4,4)

        lint    = 6

c     7-point nodal integration

      elseif(l.eq.7) then
        gg = 1.d0/15.d0
        do j = 1,3
          do i = 1,3
            el(i,j) = 0.0d0
          end do ! i
          el(j,j) = 1.00d0
          el(4,j) = 0.05d0
        end do ! j
        el(1,4) = 0.5d0
        el(2,4) = 0.5d0
        el(3,4) = 0.0d0
        el(4,4) = 2.0d0/15.0d0

        el(1,5) = 0.0d0
        el(2,5) = 0.5d0
        el(3,5) = 0.5d0
        el(4,5) = el(4,4)

        el(1,6) = 0.5d0
        el(2,6) = 0.0d0
        el(3,6) = 0.5d0
        el(4,6) = el(4,4)

        el(1,7) = one3
        el(2,7) = one3
        el(3,7) = one3
        el(4,7) = 0.45d0

        lint    = 7

      else
        if(ior.lt.0) then
          write(*,2000) l
        endif
        write(iow,2000) l
        write(ilg,2000) l
        call plstop()
      endif

c     Format

2000  format(' *ERROR* TINT2DN: Wrong quadrature, nel =',i3)

      end
