c$Id:$
      subroutine umesh8(tx,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add tx(*) to argument list                       26/09/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Dummy user input routine

c      Inputs:
c         tx(*)  - Command line input data
c         prt    - Flag, output results if true

c      Outputs:
c         none   - Users are responsible for generating outputs
c                  through common blocks, etc.  See programmer
c                  manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'umac1.h'

      logical   prt,pcomp
      character tx(*)*15

c     Set command

      if(pcomp(uct,'mes8',4)) then      ! Usual    form
c       uct = 'name'                    ! Specify 'name'
      elseif(ucount) then               ! Count elements and nodes

      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation

      endif

      end
