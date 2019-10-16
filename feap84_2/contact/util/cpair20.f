c$Id:$
      subroutine cpair20 (npair,ncdim,cp0,tydat)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor           April 10, 1996            1.0

c      Acronym: Contact PAIR description

c      Purpose: User routine for pair definitions

c      Inputs :
c         npair   - Pair number
c         ncdim   - Pair dimension
c         tydat(*)- Type data

c      Outputs:
c         cp0(*)  - Pair table
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_pair.h'

      integer   npair,ncdim
      real*8    cp0(nr0,n0c3:*),tydat(*)

      save

      call cdebug0 ('    cpair20',npair)

      write(*,3001)

3001  format (' WARNING - dummy contact pair CPAIR20 called')

      end
