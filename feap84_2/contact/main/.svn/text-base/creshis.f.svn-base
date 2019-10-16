c$Id:$
      subroutine creshis (ch1,ch2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact RESet HIStory variables

c      Purpose: Initialize t_n+1 contact history variables from final
c               value of variables at t_n

c      Inputs :
c         ch2(*)  - Contact history variables (current)

c      Outputs:
c         ch1(*)  - Contact history variables (old)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_keyh.h'

      integer   k
      real*8    ch1(*),ch2(*)

      save

      call cdebug0 ('  creshis',-1)

c     Copy history variables from vector ch2 to vector ch1

      do k = 1,tlch1
        ch1(k) = ch2(k)
      end do

      end
