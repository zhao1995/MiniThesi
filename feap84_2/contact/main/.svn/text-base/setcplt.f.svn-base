c$Id:$
      subroutine setcplt (ifplt,ifsurf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: SET Contact PLoT

c      Purpose: define pairs and surfaces that have to be plotted

c      Inputs :

c      Outputs:
c         ifplt   - PLoT flag
c         ifsurf  - SURFace flag (slave, master, both)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_pair.h'
      include  'ldata.h'
      include  'chdata.h'

      logical   ifplt,  errck,vinput,whfl,pcomp
      integer   ifsurf, firstpair,lastpair,kp
      real*8    td(3)

      save

      call cdebug0 ('        setcplt',-1)

c     Get command parameters

      if (pcomp(yyy(1:4),'pair',4)) then
        errck = vinput(yyy(16:60),45,td,3)
      else
        errck = vinput(lzz(l),50,td,3)
      endif
      firstpair = nint(td(1))
      lastpair  = nint(td(2))
      ifsurf    = nint(td(3))

c     Select pair

      if (firstpair.eq.0) firstpair = 1
      if (lastpair.eq.0)  lastpair  = firstpair

c     Check if selected pair is the current one

      whfl = .true.
      kp = firstpair - 1
      do while (whfl)
        kp = kp+1
        if (kp.gt.lastpair) then
          ifplt = .false.
          whfl  = .false.
        elseif (kp.eq.rnpair) then
          ifplt = .true.
          whfl  = .false.
        endif
      end do

      end
