c$Id:$
      subroutine setcprt (ifprt,firstel,lastel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: SET Contact PRinTout

c      Purpose: define pairs and element that have to be printed

c      Inputs :

c      Outputs:
c         ifprt   - PRinTout control flag
c         firstel - FIRST ELement to print
c         lastel  - LAST ELement to print
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_geom.h'
      include  'c_pair.h'
      include  'chdata.h'
      include  'ldata.h'

      logical   ifprt,  errck,vinput
      integer   firstel,lastel, selpair
      real*8    td(3)

      save

c     Get command parameters

      errck = vinput(lzz(l),50,td,3)

      selpair = nint(td(1))
      firstel = nint(td(2))
      lastel  = nint(td(3))

c     Select pair

      if (selpair.eq.0) selpair = 1

c     Check if selected pair is the current one

      if (selpair.eq.rnpair) then
        ifprt = .true.

c       Select elements

        if (firstel.le.0)   firstel = 1
        if (lastel.le.0)    lastel  = firstel
        if (lastel.gt.neps1) lastel = neps1
      else
        ifprt = .false.
      endif

      end
