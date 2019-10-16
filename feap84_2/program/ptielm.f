c$Id:$
      subroutine ptielm(ix,ip,nen,nen1,numel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Eliminate deleted node numbers from element lists after
c               ties

c      Inputs:
c        ix(nen1,*)  - Original connection list
c        ip(*)       - Merge node number list
c        nen         - Maximum number of nodes on an element
c        nen1        - Dimension of ix array
c        numel       - Number of elements

c      Outputs:
c        ix(nen1,*)  - Corrected connection list
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      integer    nen,nen1,numel, i,j,k
      integer    ix(nen1,numel),ip(*)

      save

c     Eliminate node j from element connections for solution

      do j = 1,numel
        do i = 1,nen
          k = abs(ix(i,j))
          if(k.gt.0) ix(i,j) = ip(k)
        end do ! i
      end do ! j

      end
