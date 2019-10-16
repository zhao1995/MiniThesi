c$Id:$
      subroutine pfacev(ix,x,ndm,iln,ct,ip,nface)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set xl(k,4) to xl(k,1) when nnode = 3            11/09/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot visible exterior faces of 3-d objects

c      Inputs:
c         ix(4)      - List of face nodes
c         x(ndm,*)   - Nodal coordinates
c         iln(2)     - Line type data
c         ct         - Also plot back faces when > 0
c         ip(8,*)    - Sort data for hidden surface representations
c         nface      - Number of faces

c      Outputs:
c         none       - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comblk.h'
      include  'pointer.h'

      logical   visbl
      integer   j, k, ndm, ip, nface, nnode
      integer   ix(4),iln(2),ilnd(2)
      real*8    x(ndm,*),xl(3,4),ct

      save

      data      ilnd(1)/1/

      ilnd(2) = iln(2)

c     Plot face

      if(ix(4).eq.0) then
        nnode = 3
      else
        nnode = 4
      endif
      do j = 1,nnode
        do k = 1,3
          xl(k,j) = x(k,ix(j))
        end do ! k
      end do ! j
      if(nnode.eq.3) then
        do k = 1,3
          xl(k,4) = xl(k,1)
        end do ! k
      endif

      if(visbl(xl) .or. (ix(1).eq.ix(4) .and. ix(2).eq.ix(3))) then

        nface = nface + 1
        do j = 1,nnode
          mr(np(66)+ix(j)-1) = 1
        end do ! j
        if(ct.gt.0.0d0) then
          call plline(iln)
          call plotl(xl(1,nnode),xl(2,nnode),xl(3,nnode),3)
          do j = 1,nnode
            call plotl(xl(1,j),xl(2,j),xl(3,j),2)
          end do ! j
        end if

      elseif(ct.gt.0.0d0) then

        ip = 0
        call plline(ilnd)
        call plotl(xl(1,nnode),xl(2,nnode),xl(3,nnode),3)
        do j = 1,nnode
          call plotl(xl(1,j),xl(2,j),xl(3,j),2)
        end do ! j

      else

        ip = 0

      endif

      end
