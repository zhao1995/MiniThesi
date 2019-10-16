c$Id:$
      subroutine opcon(numel,numcel,nen,nen1,ncen,ncen1,ix,ixc,ic,
     &                 ip,nnel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set elements connected to each node

c     Inputs:
c       numel        - Number elements
c       numcel       - Number contact elements
c       nen          - Number nodes (max) connected to element
c       nen1         - Dimension for 'ix' array
c       ncen         - Number nodes (max) connected to contact element
c       ncen1        - Dimension for 'ixc' array
c       ix(nen1,*)   - Element nodal connections
c       ixc(ncen1,*) - Element nodal connections
c       ic(*)        - Pointer for storage

c     Outputs:
c       ip(*)        - Elements connected to each node
c       nnel(*)      - Number nodes on each element
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    numel,numcel,nen,nen1,ncen,ncen1, i,ii,jj,n
      integer    ix(nen1,*),ixc(ncen1,*),ic(0:*),ip(*),nnel(*)

      save

c     Do normal elements

      do n = 1,numel
        do i = 1,nen
          if(ix(i,n).gt.0) then
            ii      = ix(i,n)
            nnel(n) = i
            do jj = ic(ii-1)+1,ic(ii)
              if(ip(jj).eq.0) then
                ip(jj) = n
                go to 100
              endif
            end do ! jj
            write(*,*) 'ERROR in POPCON: IX'
          endif
100       continue
        end do ! i
      end do ! n

c     Do contact elements

      do n = 1,numcel
        do i = 1,ncen
          if(ixc(i,n).gt.0) then
            ii            = ixc(i,n)
            nnel(n+numel) = i
            do jj = ic(ii-1)+1,ic(ii)
              if(ip(jj).eq.0) then
                ip(jj) = n + numel
                go to 200
              endif
            end do ! jj
            write(*,*) 'ERROR in POPCON: IXC'
          endif
200       continue
        end do ! i
      end do ! n

      end
