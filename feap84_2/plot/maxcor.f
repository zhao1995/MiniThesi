c$Id:$
      subroutine maxcor(x,ndm,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine max/min coordinates for a perspective
c               view selection

c      Inputs:
c         x(ndm,*)  - Nodal coordinates (may be in deformed state)
c         ndm       - Spatial dimension
c         numnp     - Number of nodes in mesh

c      Outputs:
c         none      - Output through common /pdata0/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata0.h'
      include  'pdatas.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   ndm, numnp, i,ii,n
      real*8    x(ndm,*), vm

      save

c     Determine max/min coordinates for a perspective view selection

      ii    = min(ndm,3)
      fp(1) = npty - 1
      do i = 1,ii
        vmin(i) = x(i,1)
        do n = 1,numnp
            vmin(i) = max(vmin(i),x(i,n))
        end do ! n
        vmax(i) = vmin(i)
        do n = 1,numnp
          if(mr(fp(1)+n).ge.0) then
            vmin(i) = min(vmin(i),x(i,n))
          endif
        end do ! n
      end do ! i

c     Allow for symmetries in selecting view

      do i = 1,ii
        if(isymm(i,1).ne.0) then
          vm      = min( vmin(i), -vmax(i), - vmin(i))
          vmax(i) = max( vmax(i), -vmax(i), - vmin(i))
          vmin(i) = vm
        endif
      end do ! i

      end
