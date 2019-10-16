c$Id:$
      subroutine int2dc(l,lint,sw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Integration on circular beam sections

c     Inputs:
c        l       - Method type
c                  1   : 4-point method (all interior)
c                  2   : 5-point method (center + perimeter)
c                  3   : 9-point method (center, interior, perimeter)
c                  4   :17-point method (center, interior, perimeter)

c     Outputs:

c        lint    - Number of quadrature points
c        sw(3,*) - Points and weights
c                  1,2 : Points
c                  3   : Weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'iofile.h'
      include 'pconstant.h'

      integer  l,lint
      real*8   sw(3,*), rt2

      save

c     4-Point method: R = O(h^4)

      if(l.eq.1) then

        lint = 4

        sw(1,1) =  0.5d0
        sw(2,1) =  0.5d0
        sw(3,1) =  0.25d0

        sw(1,2) =  0.5d0
        sw(2,2) = -0.5d0
        sw(3,2) =  0.25d0

        sw(1,3) = -0.5d0
        sw(2,3) = -0.5d0
        sw(3,3) =  0.25d0

        sw(1,4) = -0.5d0
        sw(2,4) =  0.5d0
        sw(3,4) =  0.25d0

c     5-Point method: R = O(h^4)

      elseif(l.eq.2) then

        lint = 5

        sw(1,1) =  1.0d0
        sw(2,1) =  0.0d0
        sw(3,1) =  0.125d0

        sw(1,2) =  0.0d0
        sw(2,2) =  1.0d0
        sw(3,2) =  0.125d0

        sw(1,3) = -1.0d0
        sw(2,3) =  0.0d0
        sw(3,3) =  0.125d0

        sw(1,4) =  0.0d0
        sw(2,4) = -1.0d0
        sw(3,4) =  0.125d0

        sw(1,5) =  0.0d0
        sw(2,5) =  0.0d0
        sw(3,5) =  0.5d0

c     9-Point method: R = O(h^6)

      elseif(l.eq.3) then

        lint = 9

        sw(1,1) =  1.0d0
        sw(2,1) =  0.0d0
        sw(3,1) =  1.d0
        sw(3,1) =  sw(3,1)/24.0d0

        sw(1,2) =  0.0d0
        sw(2,2) =  1.0d0
        sw(3,2) =  sw(3,1)

        sw(1,3) = -1.0d0
        sw(2,3) =  0.0d0
        sw(3,3) =  sw(3,1)

        sw(1,4) =  0.0d0
        sw(2,4) = -1.0d0
        sw(3,4) =  sw(3,1)

        sw(1,5) =  0.5d0
        sw(2,5) =  0.5d0
        sw(3,5) =  one6

        sw(1,6) =  0.5d0
        sw(2,6) = -0.5d0
        sw(3,6) =  sw(3,5)

        sw(1,7) = -0.5d0
        sw(2,7) = -0.5d0
        sw(3,7) =  sw(3,5)

        sw(1,8) = -0.5d0
        sw(2,8) =  0.5d0
        sw(3,8) =  sw(3,5)

        sw(1,9) =  0.0d0
        sw(2,9) =  0.0d0
        sw(3,9) =  sw(3,5)

c     17-Point method: R = O(h^6)

      elseif(l.eq.4) then

        lint = 17

c       Perimeter

        rt2  = 1.d0/sqrt2

        sw(1,1) =  1.0d0
        sw(2,1) =  0.0d0
        sw(3,1) =  1.d0
        sw(3,1) =  sw(3,1)/48.0d0

        sw(1,2) =  rt2
        sw(2,2) =  rt2
        sw(3,2) =  sw(3,1)

        sw(1,3) =  0.0d0
        sw(2,3) =  1.0d0
        sw(3,3) =  sw(3,1)

        sw(1,4) = -rt2
        sw(2,4) =  rt2
        sw(3,4) =  sw(3,1)

        sw(1,5) = -1.0d0
        sw(2,5) =  0.0d0
        sw(3,5) =  sw(3,1)

        sw(1,6) = -rt2
        sw(2,6) = -rt2
        sw(3,6) =  sw(3,1)

        sw(1,7) =  0.0d0
        sw(2,7) = -1.0d0
        sw(3,7) =  sw(3,1)

        sw(1,8) =  rt2
        sw(2,8) = -rt2
        sw(3,8) =  sw(3,1)

c       Interior Points

        sw(1, 9) =  rt2
        sw(2, 9) =  0.0d0
        sw(3, 9) =  1.d0
        sw(3, 9) =  sw(3,9)/12.d0

        sw(1,10) =  0.5d0
        sw(2,10) =  0.5d0
        sw(3,10) =  sw(3,9)

        sw(1,11) =  0.0d0
        sw(2,11) =  rt2
        sw(3,11) =  sw(3,9)

        sw(1,12) = -0.5d0
        sw(2,12) =  0.5d0
        sw(3,12) =  sw(3,9)

        sw(1,13) = -rt2
        sw(2,13) =  0.0d0
        sw(3,13) =  sw(3,9)

        sw(1,14) = -0.5d0
        sw(2,14) = -0.5d0
        sw(3,14) =  sw(3,9)

        sw(1,15) =  0.0d0
        sw(2,15) = -rt2
        sw(3,15) =  sw(3,9)

        sw(1,16) =  0.5d0
        sw(2,16) = -0.5d0
        sw(3,16) =  sw(3,9)

c       Center

        sw(1,17) =  0.0d0
        sw(2,17) =  0.0d0
        sw(3,17) =  one6

      else

        write(ilg,2000) l
        write(iow,2000) l
        if(ior.lt.0) then
          write(*,2000) l
        endif
        call plstop()

      endif

2000  format(' *ERROR* INT2DC: Circular Integration Type',i3,
     &       ' not Impemented')

      end
