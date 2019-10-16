c$Id:$
      subroutine crestrt (csw,ch1,ch2,ch3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0
c               Robert L. Taylor         January 19, 2003            1.1

c      Acronym: Contact RESTaRT

c      Purpose: Read or write data for restart (No checks performed)

c      Inputs :

c      Outputs:
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_keyh.h'
      include  'iodata.h'

      integer   csw
      real*8    ch1(tlch1),ch2(tlch1),ch3(tlch3)

      save

      call cdebug0 ('  crestrt',-1)

c     Read data

      if (csw.eq.306) then

        read (ios) ch1,ch2,ch3

c     Write data

      elseif (csw.eq.307) then

        write (ios) ch1,ch2,ch3

      endif

      end
