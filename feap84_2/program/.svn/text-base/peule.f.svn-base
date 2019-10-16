c$Id:$
      subroutine peule(ix,nen,leule,euler,nrot)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set value of Euler angle for each element node

c      Inputs:
c         ix(nen)    - Element nodal connection list
c         nen        - Number of nodes connected to element
c         euler(3,*) - Nodal Euler angle array

c      Outputs:
c         leule(3,*) - Element Euler angle array
c         nrot       - Number of nodes needing modification
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nen,nrot,ii,n,j
      integer   ix(nen)
      real*8    leule(3,nen),euler(3,*)

      save

c     Set up table of Euler angles

      nrot = 0
      do n = 1,nen
        do j = 1,3
          leule(j,n) = 0.0d0
        end do ! j
        ii = ix(n)
        if (ii.gt.0) then
          if (euler(1,ii).ne.0.0d0 .or. euler(2,ii).ne.0.0d0
     &                             .or. euler(3,ii).ne.0.0d0) then
            do j = 1,3
              leule(j,n) = euler(j,ii)
            end do ! j
            nrot = nrot + 1
          endif
        endif
      end do ! n

      end
