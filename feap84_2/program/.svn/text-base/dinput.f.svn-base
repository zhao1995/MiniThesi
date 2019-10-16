c$Id:$
      subroutine dinput(d,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c      N.B.  This routine is superceded by 'pinput'

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Data input subprogram for real values

c      Inputs:
c         nn    - Number of data items to input

c      Outputs:
c         d(nn) - Values for nn items input
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'errchk.h'

      logical   pinput
      integer   nn
      real*8    d(nn)

      save

c     Input nn items

      errck = pinput(d,nn)

      end
