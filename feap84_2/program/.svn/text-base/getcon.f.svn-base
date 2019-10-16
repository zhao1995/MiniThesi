c$Id:$
      subroutine getcon(epmac)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform installation parameter computations

c      Inputs:

c      Outputs:
c         epmac   - Smallest number that can be added to 1.0
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    epmac

      save

c     Compute machine epsilon estimate

      epmac = 1.0d0
100   epmac = epmac*0.5d0
      if(1.d0 .ne. 1.d0 + epmac) go to 100
      epmac = 2.0d0*epmac

      end
