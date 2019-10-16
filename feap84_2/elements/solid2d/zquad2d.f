c$Id:$
      subroutine zquad2d(lint, linr, sg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Enhanced mode shaped functions

c      Inputs:
c        lint    - Number of quadrature points

c      Outputs:
c        linr    - Number of projection points
c        sg(2,*) - Location of projection points
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   lint, linr
      real*8    sg(2,*)

      save

c     Linear elements

      if(lint.eq.4) then

        sg(1,1) = 0.0d0     
        sg(2,1) = 0.0d0     

        linr    = 1

c     Quadratic elements

      elseif(lint.eq.8 .or. lint.eq.9) then
        sg(1,1) = -1.d0/sqrt(1.8d0)
        sg(2,1) =  sg(1,1)

        sg(1,2) = -sg(1,1)
        sg(2,2) =  sg(2,1)

        sg(1,3) = -sg(1,1)
        sg(2,3) = -sg(2,1)

        sg(1,4) =  sg(1,1)
        sg(2,4) = -sg(2,1)

        linr    = 4

c     Cubic elements

      elseif(lint.eq.12 .or. lint.eq.16) then

        linr    = 9

      endif

      end
