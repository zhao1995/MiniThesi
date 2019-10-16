c$Id:$
      subroutine optic(numnp,numel,numcel,nen,nen1,ncen,ncen1,ix,ixc,
     &                 ic,kp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute number elements attached to each node

c     Inputs:
c       numnp      - Number nodes
c       numel      - Number elements
c       numcel     - Number contact elements
c       numel      - Number material sets
c       nen        - Number nodes (max) connected to element
c       nen1       - Dimension for 'ix' array
c       ncen       - Number nodes (max) connected to contact element
c       ncen1      - Dimension for 'ixc' array
c       ix(nen1,*) - Element nodal connections

c     Outputs:
c       ic(*)      - Pointer for storage
c       kp         - Storage required for 'ip' array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    numnp,numel,numcel,nen,nen1,ncen,ncen1, n,i,kp
      integer    ix(nen1,*),ixc(ncen1,*), ic(0:numnp)

      save

c     Initialize pointer array

      do n = 0,numnp
        ic(n) = 0
      end do ! n

c     Count occurances of nodes in elements

      do n = 1,numel
        do i = 1,nen
          if(ix(i,n).gt.0) then
            ic(ix(i,n)) = ic(ix(i,n)) + 1
          end if
        end do ! i
      end do ! n

      do n = 1,numcel
        do i = 1,ncen
          if(ixc(i,n).gt.0) then
            ic(ixc(i,n)) = ic(ixc(i,n)) + 1
          end if
        end do ! i
      end do ! n

c     Convert to pointers

      do n = 1,numnp
        ic(n) = ic(n) + ic(n-1)
      end do ! n

c     Maximum array side needed

      kp = ic(numnp)

      end
