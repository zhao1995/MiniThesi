c$Id:$
      subroutine inchv1 (ch1,ch2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: INit Contact HIstory Variables

c      Purpose: set data at time t=0 which certainly formulations need

c      Inputs :
c         ch2(*)  - Contact history variables (current)

c      Outputs:
c         ch1(*)  - Contact history variables (old)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_keyh.h'
      include  'c_pair.h'

      real*8    ch1(*),ch2(*), masts,istgn,d21,csi,gt

      save

      call cdebug0 ('      inchv1',-1)

c     variables for Coulomb friction

      if (iffric.eq.1) then
        masts = ch2(p1(1))
        istgn = ch2(p1(4))
        d21   = ch2(p1(7))
        csi   = ch2(p1(8))
        gt    = ch2(p1(10))

        ch1(p1(1))  = masts
        ch1(p1(4))  = istgn
        ch1(p1(7))  = d21
        ch1(p1(8))  = csi
        ch1(p1(10)) = gt
        ch1(p1(17)) = csi
        ch1(p1(58)) = -2
      endif

      if (ifsolm.eq.2) then
        ch1(p1(21)) = 0.0d0
      endif

      end
