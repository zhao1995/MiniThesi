c$Id:$
      subroutine pltmv(pl,ipl,u,nplts,save)


c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Pick component from second entry, not first      24/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Extract information for time history from solution
c               array

c      Inputs:
c         ipl(2,*) - Location of information in u array
c         u(*)     - Array of values
c         nplts    - Number of items to extract
c         save     - Sign of quantity to save

c      Outputs:
c         pl(*)    - Array of time history information
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nplts,j,n, ipl(2,nplts)
      real*8    save, pl(nplts),u(*)

      save

      do n = 1,nplts
        j = ipl(2,n)
        if(j.gt.0) then
          pl(n) = u(j)*save
        end if
      end do ! n

      end
