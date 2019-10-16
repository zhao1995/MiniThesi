c$Id:$
      subroutine int1d(l,sw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Start 5 point at negative values                 08/09/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Gauss quadrature for 1-d element

c      Inputs:
c         l     - Number of points

c      Outputs:
c         sw(1,*) - Gauss points
c         sw(2,*) - Gauss weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   l
      real*8    sw(2,*), t

      save

      if(l.eq.1) then

        sw(1,1) = 0.0d0
        sw(2,1) = 2.0d0

      elseif(l.eq.2) then

        sw(1,1) = -sqt13                     ! sqrt(1/3)
        sw(1,2) = -sw(1,1)
        sw(2,1) = 1.0d0
        sw(2,2) = 1.0d0

      elseif(l.eq.3) then

        sw(1,1) = -sqtp6                    ! sqrt(0.6)
        sw(1,2) = 0.0d0
        sw(1,3) = -sw(1,1)
        sw(2,1) = five9
        sw(2,2) = eight9
        sw(2,3) = sw(2,1)

      elseif(l.eq.4) then

        t       =  sqt48                    ! sqrt(4.8)
        sw(1,1) = -sqrt((3.d0+t)/7.d0)
        sw(1,2) = -sqrt((3.d0-t)/7.d0)
        sw(1,3) = -sw(1,2)
        sw(1,4) = -sw(1,1)
        t       =  one3/t                   ! 1/3/t
        sw(2,1) =  0.5d0 - t
        sw(2,2) =  0.5d0 + t
        sw(2,3) =  sw(2,2)
        sw(2,4) =  sw(2,1)

      elseif(l.eq.5) then


        t       =  1120.d0
        t       =  sqrt(t)

        sw(1,1) = -sqrt((70.d0 + t)/126.d0)
        sw(1,2) = -sqrt((70.d0 - t)/126.d0)
        sw(1,3) =  0.0d0
        sw(1,4) = -sw(1,2)
        sw(1,5) = -sw(1,1)

        sw(2,1) =  (21.d0*t + 117.6d0)/(t*(70.d0 + t))
        sw(2,2) =  (21.d0*t - 117.6d0)/(t*(70.d0 - t))
        sw(2,3) =  2.d0*(1.d0 - sw(2,1) - sw(2,2))
        sw(2,4) =  sw(2,2)
        sw(2,5) =  sw(2,1)

c     Compute points and weights

      else

        call gausspw(l,sw)

      endif

      end
