c$Id:$
      subroutine cextvar (ch1,ch2,ch3,chvec,chvar,dr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise            July 10, 1996            1.0

c      Acronym: Contact EXTract history VARiable

c      Purpose: Extract history variable for pair and store on file

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_keyh.h'

      integer   chvec,chvar, kv
      real*8    ch1(lh1,*),ch2(lh1,*),ch3(lh3,*),dr(*)

      save

      if (chvec.eq.1) then
        do kv = 1, nset
          dr(kv) = ch1(p1(chvar),kv)
        end do
      elseif (chvec.eq.2) then
        do kv = 1, nset
          dr(kv) = ch2(p1(chvar),kv)
        end do
      else
        do kv = 1, nset
          dr(kv) = ch3(p3(chvar),kv)
        end do
      endif

      end
