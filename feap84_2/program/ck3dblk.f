c$Id:$
      subroutine ck3dblk (ix, xl, nel, ndm, shp, err )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check 3d block for input sequence errors
c               Write message to file on errors located.

c      Inputs:
c         ix(*)     - List of nodes for block
c         xl(ndm,*) - Coordinate array
c         nel       - Number of nodes on block
c         ndm       - Spatial dimension of mesh

c      Outputs:
c         None

c      Scratch:
c         shp(*)    - Array to store shape functions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   err
      integer   ndm, nel, i, l, ineg, ix(*), ic(8)
      real*8    detj, rst(3,8), xl(ndm,*), shp(*), jac(8)

      save

      data      rst/-1.d0,-1.d0,-1.d0,   1.d0,-1.d0,-1.d0,
     &               1.d0, 1.d0,-1.d0,  -1.d0, 1.d0,-1.d0,
     &              -1.d0,-1.d0, 1.d0,   1.d0,-1.d0, 1.d0,
     &               1.d0, 1.d0, 1.d0,  -1.d0, 1.d0, 1.d0/

c     Compute jacobian at each corner of element

      if(nel.ge.8) then
        ineg = 0
        err  = .false.
        do l = 1,min(8,nel)
          call bjac3d ( rst(1,l) , xl, ndm, shp, detj )
          if(detj.le.0.0d0) then
            ineg      = ineg + 1
            ic(ineg)  = ix(l)
            jac(ineg) = detj
          endif
        end do ! l
        if(ineg.gt.0) then
          err = .true.
          write(iow,2001) (ic(i),jac(i),i=1,ineg)
          if(ior.lt.0) write(*,2001) (ic(i),jac(i),i=1,ineg)
          call iprint(ix,1,8,1,'IXL-nodes')
          call mprint(xl,3,8,ndm,'XL-coord')
        endif
        if(ineg.eq.nel) then
          write(iow,2002)
          if(ior.lt.0) then
            write(*,2002)
          endif
        endif
      endif

2001  format(/5x,'*ERROR*: BLOCk has zero or negative jacobian':,
     &        ' at block nodes:'/
     &      (10x,'Node =',i3,' Jacobian =',1p,1e12.5))

2002  format(/5x,'To correct reverse node sequence on block')

      end
