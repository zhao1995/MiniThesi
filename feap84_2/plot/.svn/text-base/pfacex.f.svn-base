c$Id:$
      subroutine pfacex(il,ix,ixf,nen,nen1,nf,n)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set up nodes for faces

c      Inputs:
c         il(*)   - Location of face nodes on element
c         ix(*)   - Node numbers on elements
c         nen     - Number of nodes on element
c         nen1    - Location of material set number on element
c         n       - Element number

c      Outputs:
c         ixf(*)  - Face nodes
c         nf      - Face number
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nen,nen1,nf,n, j
      integer   il(4),ix(nen1),ixf(7)

      save

c     Set face nodes

      do j = 1,min(4,nen)
        ixf(j) = ix(il(j))
      end do ! j

c     Set region and material number

      ixf(5) = n
      ixf(6) = ix(nen1-1)
      ixf(7) = ix(nen1)
      nf = nf + 1

      end
