c$Id:$
      subroutine chkblk(y,n0,nt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Add zero entry to character array which has blank
c               field.

c      Inputs:
c         y(*) - array to check
c         n0   - Field width of data
c         nt   - Size of array to check

c      Outputs:
c         y(*) - Blank fields have zero (0) added to field 'n0'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   pcomp
      integer   n0,nt,n
      character y*(*)

      save

c     Patch for plot device type if necessary

      if(pcomp(y,'tekt',4)) then
        y(16:19) = y(27:30)
      endif

c     Add character if y(nt) is blank

      do n = n0,nt,n0
        if(y(n:n).eq.' ') y(n:n) = '0'
      end do ! n

      end
