c$Id:$
      subroutine cpoutm (cs0,cm0,cp0,ics,hic)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    26/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Control for output of parallel contact data to file
c              Dummy in serial version

c     Inputs:
c       cs0(*)    - Surface  table
c       cm0(*)    - Material table
c       cp0(*)    - Pair     table
c       ics(*)    - Contact facet data
c       hic(*)    - History correspondence table

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    ics(*), hic(*)
      real*8     cs0(*), cm0(*), cp0(*)

      end




