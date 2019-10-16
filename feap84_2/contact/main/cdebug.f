c$Id:$
      subroutine cdebug (nn,vect,strin)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact DEBUG

c      Purpose: dump requested storage area

c      Inputs :
c         nn      - Set to be dumped (0 = all)
c         vect    - VECTor to be dumped
c         strin   - identification STRINg for call

c      Outputs:
c                 - In file 'Cdebug'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_ccp.h'
      include  'c_chp.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'print.h'
      include  'comblk.h'

      logical   screen
      character vect*(*),strin*(*)
      integer   nn, debf

      save

c     Set debug file for memory dump

      if (vect.eq.'debf') then
        if (nn.ne.0) then
          debf = nn
        else
          debf = 98
        endif
        open (unit = debf, file = 'Cdebug', status = 'unknown')

c     Perform dump only if debug mode is active
c     or request from a 'show' command
c     WARNING dump in file is checked with ifdb in interactive mode

      elseif ((ifdb) .or. (strin.eq.'interactive')) then

c       Set switch for screen output

        if (strin.eq.'interactive') then
          screen = .true.
        else
          screen = .false.
        endif

c       Print label

        if (ifdb) then
          write (debf,10) strin
        endif
        if (screen) then
          write (iow,10) strin
        endif
        if ((ior.lt.0) .and. (screen)) then
          write (*,10) strin
        endif

c       Perform dump

        call dumpmem (debf,nn,vect,hr(ccp(1)),mr(np(133)),hr(ccp(2)),
     &                hr(np(132)),hr(ccp(3)),hr(chp(1)),hr(chp(2)),
     &                hr(chp(3)),c_lp1,c_lp3,mr(np(134)),screen)
      endif

c     Formats

10    format (//a60)

      end
