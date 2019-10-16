c$Id:$
      logical function pinput(d,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input routine for real data.: returns true on error

c      Inputs:
c         nn     - Number of data items to extract from input string
c                  N.B. Input performed by this function

c      Outputs:
c         d(*)   - Values of data input
c         pinput - Flag, returns true if error occurs during input
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character tl*15
      logical   tinput
      integer   nn
      real*8    d(nn)

      save

c     Use routine TINPUT to perform the inputs

      pinput = tinput(tl,0,d,nn)

      end
