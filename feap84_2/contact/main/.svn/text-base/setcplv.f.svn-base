c$Id:$
      subroutine setcplv (ifplt,chvec,chvar)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: SET Contact PLot Variable

c      Purpose: Define pairs & history variable that has to be plotted

c      Inputs :

c      Outputs:
c         ifplt   - PLoT flag
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_pair.h'
      include  'ldata.h'
      include  'chdata.h'

      logical   ifplt,  errck,vinput,pcomp
      integer   chvec,chvar, selpair
      real*8    td(3)

      save

      call cdebug0 ('        setcplt',-1)

c     Get command parameters

      if (pcomp(yyy(1:4),'cvar',4)) then
        errck = vinput(yyy(16:80),50,td,3)
      else
        errck = vinput(lzz(l),50,td,3)
      endif

c     Set contact variable number

      chvar = max(1,nint(td(1)))

c     Set contact vector number: 1 = ch1; 2 = ch2; 3 = ch3

      chvec = min(3,nint(td(2)))
      if (chvec.eq.0) then
        chvec = 2
      endif

c     Select pair

      selpair = nint(td(3))
      if (selpair.eq.0) then
        ifplt = .true.
      elseif (selpair.eq.rnpair) then
        ifplt = .true.
      else
        ifplt = .false.
      endif

      end
