c$Id:$
      subroutine ploa1(time,dt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Calculate load level for arc length method

c      Inputs:
c         time     - Current solution time
c         dt       - Current solution time increment

c      Outputs:
c         rlnew    - Arc length load level for step
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'arclel.h'
      include  'arcler.h'

      real*8    time,dt

      save

c     For restart only

      if (refl) then
        refl =.false.
      endif

c     Check if new time step

      if(arcf .and. time .ne. timold) then

c       Set new load level

        rlnew = rlnew + dt

      elseif(.not. arcf) then

        rlnew = 1.d0

      endif

      end
