c$Id:$
      subroutine autdt( dtold, dtnew, upflag )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Adjust time step using logarithmic increment adjustment

c      Inputs:
c         dtold  - size of current time step
c         upflag - .true. Use upfact in adjusting time step
c                  .false. Do not use upfact in adjustment

c      Outputs:
c         dtnew  - size of adjusted time step
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'auto1.h'
      include  'iofile.h'
      include  'print.h'

      logical   upflag
      real*8    dtold, dtnew, upfact, dnfact

      save

      data      upfact, dnfact/ 0.20d0, 0.20d0/

c     Set new time step value

      if(upflag) then

        dtnew = min(10.d0**(log10(dtold) + upfact),dtmax)

      else

        dtnew = max(10.d0**(log10(dtold) - dnfact),dtmin)

      endif

      if(dtold.ne.dtnew .and. prnt) then
        write(iow,*) 'dt-old =',dtold
        write(iow,*) 'dt-new =',dtnew
        if(ior.lt.0) then
          write(*,*) 'dt-old =',dtold
          write(*,*) 'dt-new =',dtnew
        endif
      endif

      end
