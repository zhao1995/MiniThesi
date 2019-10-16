c$Id:$
      subroutine skipcont ()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: SKIP CONTact data

c      Purpose: Skip input data if contact os "off"

c      Inputs:
c         From read of data file(s)

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   errck,gettxtd,pcomp
      character cc*4,tx(2)*15
      real*8    td(15)

      save

      call cdebug0 ('  skipcont',-1)

c     Scan and skip input data

      errck = gettxtd(tx,2,td,0,'skip')
      cc = tx(1)(1:4)

      do while (.not.pcomp(cc,'end',3))
        errck = gettxtd(tx,2,td,0,'skip')
        cc    = tx(1)(1:4)
      end do

      end
