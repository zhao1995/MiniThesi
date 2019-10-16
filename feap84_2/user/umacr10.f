c$Id:$
      subroutine umacr10(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove 'prt' from argument list                  09/07/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  User interface for adding solution command language
c                instructions.

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters

c      Outputs:
c         N.B.  Users are responsible for command actions.  See
c               programmers manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'umac1.h'

      logical   pcomp
      character lct*15
      real*8    ctl(3)

      save

c     Set command word

      if(pcomp(uct,'ma10',4)) then      ! Usual    form
c       uct = 'name'                    ! Specify 'name'
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation

      endif

      end
