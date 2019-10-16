c$Id:$
      function chp(nchv)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronim: Contact History Pointer

c      Purpose: Set the appropriate pointer for the history vectors

c      Inputs :
c         nchv    - # of the history vector

c      Outputs:
c         chp     - absolute pointer
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_keyh.h'
      include  'c_chp.h'
      include  'iofile.h'
      include  'pointer.h'

      integer   nchv

      save

c     Pointer for the contact history vector CH1

      if (nchv.eq.1) then
        chp = np(135)

c     Pointer for the contact history vector CH2

      elseif (nchv.eq.2) then
        chp = np(135) + tlch1

c     Pointer for the contact history vector CH3

      elseif (nchv.eq.3) then
        chp = np(135) + 2*tlch1

      else
        chp = 1
        write(ilg,3000) nchv
        call plstop()
      endif

c     Format

3000  format(' *ERROR* in CHP:  Attempt to assign pointer type =',i4)

      end
