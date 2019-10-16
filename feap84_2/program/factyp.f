c$Id:$
      integer function factyp(ix, nec)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1.  Add option for tetrahedral elements             02/09/2007
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Set face type for element interfaces
c     Inputs:
c       ix(*)  -  Element connection list

c     Outputs:
c       factyp - Type of interface to check
c       nec    - Number of interface to check
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit    none

      include    'cdata.h'
      include    'sdata.h'

      integer     ix(*),i,nec

      save

      do i = nen,2,-1
        if(ix(i).ne.0) exit
      end do ! i

      if(i.le.2 .or. ndm.eq.1) then
        factyp = 1
        nec    = min(i,2)
      elseif(i.le.3 .or. (i.eq.4 .and. ndm.eq.2)) then
        factyp = 2
        nec    = i
      elseif((i.eq. 4 .or. i.eq.10 .or. i.eq.11 .or.
     &        i.eq.14 .or. i.eq.15) .and. ndm.eq.3) then
        factyp = 5
        nec    = 4
      elseif(ndm.eq.2) then
        factyp = 2
        nec    = 4
      else
        factyp = 3
        nec    = 6
      endif

      end
