c$Id:$
      subroutine ckbrk8 ( n, ix, xl, ndm, nel, shp )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check 8-node brick for bad data.
c               Write message to file on errors located.

c      Inputs:
c         n         - Number of element being checked
c         ix(*)     - List of nodes for element
c         xl(ndm,*) - Coordinate array
c         ndm       - Spatial dimension of mesh
c         nel       - Number of nodes on element

c      Outputs:
c         None

c      Scratch:
c         shp(*)    - Array to store shape functions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'fdata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   ndm, nel, i, l, n, ineg, ix(*), ic(16)
      real*8    detj, rst(3,8), xl(ndm,*), shp(*)

      save

      data      rst/-1.d0,-1.d0,-1.d0,   1.d0,-1.d0,-1.d0,
     &               1.d0, 1.d0,-1.d0,  -1.d0, 1.d0,-1.d0,
     &              -1.d0,-1.d0, 1.d0,   1.d0,-1.d0, 1.d0,
     &               1.d0, 1.d0, 1.d0,  -1.d0, 1.d0, 1.d0/

c     Check element for input errors

      ineg = 0
      do l = 1,min(8,nel)
        if(ix(l).gt.0) then
          if(mr(npty+ix(l)-1).lt.0) then
            ic(ineg+1) = l
            ic(ineg+2) = abs(ix(l))
            ineg = ineg + 2
          endif
        endif
      end do ! l

c     Node numbering errors

      if(ineg.gt.0) then
        write(iow,2000) n,(ic(i),i=1,ineg)
        if(ior.lt.0) write(*,2000) n,(ic(i),i=1,ineg)

c     Compute jacobian at each corner of element

      else
        do l = 1,min(8,nel)
          call bjac3d ( rst(1,l) , xl, ndm, shp, detj )
          if(detj.le.0.0d0) then
            ic(ineg+1) = l
            ic(ineg+2) = abs(ix(l))
            ineg = ineg + 2
          endif
        end do ! l
        if(ineg.gt.0 .and. pfr) then
          write(iow,2001) n,(ic(i),i=1,ineg)
          if(ior.lt.0) write(*,2001) n,(ic(i),i=1,ineg)
        endif
        if(nel.eq.8 .and. ineg.eq.2*nel) then
          do l = 1,4
            ineg    = ix(l)
            ix(l  ) = ix(l+4)
            ix(l+4) = ineg
          end do ! l
        end if
      endif

2000  format(' >Element',i9,' coordinates not input for nodes:'/
     &      ('                Local =',i3,' Global =',i9))

2001  format(' >Element',i9,' has negative jacobian at nodes:'/
     &      ('                Local =',i3,' Global =',i9))

      end
