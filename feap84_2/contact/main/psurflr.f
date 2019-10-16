c$Id:$
      subroutine psurflr(ics,neps,nepsr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Check for closed surfaces and reset number facets

c     Inputs:
c       ics(*)     - List of original nodes for each surface.

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  neps,nepsr
      integer  ics(neps)

c     Check if ends have same node number

      if(ics(1).eq.ics(neps)) then
        nepsr = neps - 1
      else
        nepsr = neps
      endif

      end
