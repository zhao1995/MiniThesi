c$Id:$
      subroutine pltriad(ix,nen,ltriad,triad,nrot)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set value of Euler angle for each element node

c      Inputs:
c         ix(nen)      - Element nodal connection list
c         nen          - Number of nodes connected to element
c         triad(3,3,*) - Nodal triad array

c      Outputs:
c         ltriad(3,3,*) - Elemet triad array
c         nrot          - Number of nodes needing modification
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nen,nrot,ii,n
      integer   ix(nen)
      real*8    ltriad(3,3,nen),triad(3,3,*)

      save

c     Set up table of local triad arrays

      nrot = 0
      do n = 1,nen
        ltriad(:,:,n) = 0.0d0
        ii = ix(n)
        if (ii.gt.0) then
          if(triad(1,1,ii).ne.-100.d0) then
            nrot = nrot + 1
          endif
          ltriad(:,:,n) = triad(:,:,ii)
        endif
      end do ! n

      end
