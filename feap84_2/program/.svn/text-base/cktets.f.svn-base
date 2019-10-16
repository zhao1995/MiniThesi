c$Id:$
      subroutine cktets ( n, ix, xl, ndm, nel, shp )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check tetrahedra for valid data specifications.

c      Inputs:
c         n         - Element number being checked
c         ix(*)     - Number of nodes connected to element
c         xl(ndm,*) - Array of element nodal coordinates
c         ndm       - Spatial dimension of mesh
c         nel       - Number of nodes on this element

c      Outputs:
c         None

c      Scratch:
c         shp(*)    - Shape function array storage
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   ndm, nel, i, l, n, ineg, ix(*), ic(2,11)
c     integer   ll
      real*8    detj, rst(3,11), xl(ndm,*), shp(*), jac(11), jacm
c     real*8    sg(5,16)

      save

      data      rst/ 1.00d0, 0.00d0, 0.00d0,   0.00d0, 0.00d0, 1.00d0,
     &               0.00d0, 1.00d0, 0.00d0,   0.00d0, 0.00d0, 0.00d0,
     &               0.50d0, 0.00d0, 0.50d0,   0.00d0, 0.50d0, 0.50d0,
     &               0.50d0, 0.50d0, 0.00d0,   0.50d0, 0.00d0, 0.00d0,
     &               0.00d0, 0.00d0, 0.50d0,   0.00d0, 0.50d0, 0.00d0,
     &               0.25d0, 0.25d0, 0.25d0/

c     Check element for input errors

      ineg = 0
      do l = 1,nel
        if(ix(l).gt.0) then
          if(mr(np(190)+ix(l)-1).lt.0) then
            ineg       = ineg + 1
            ic(1,ineg) = l
            ic(2,ineg) = abs(ix(l))
          endif
        endif
      end do ! l

c     Node numbering errors

      if(ineg.gt.0) then
        write(iow,2000) n,(ic(1,i),ic(2,i),i=1,ineg)
        if(ior.lt.0) then
          write(*,2000) n,(ic(1,i),ic(2,i),i=1,ineg)
        endif

c     Compute jacobian at each corner of element

      else
c       call tint3d(4,ll,sg)
        ineg = 0
        jacm = 0.0d0
        do l = 1,nel
c         call tjac3d (  sg(1,l) , xl, ndm, nel, shp, detj )
          call tjac3d ( rst(1,l) , xl, ndm, nel, shp, detj )
          jacm = max(jacm,abs(detj))
          if(detj.le.0.0d0) then
            ineg       = ineg + 1
            ic(1,ineg) = l
            ic(2,ineg) = abs(ix(l))
            jac(ineg)  = detj
          endif
        end do ! l
        if(ineg.gt.0) then
          write(iow,2001) n,(ic(1,i),ic(2,i),jac(i),i=1,ineg)
          write(iow,*   ) '        Max Det J =',jacm
          if(ior.lt.0) then
            write(*,2001) n,(ic(1,i),ic(2,i),jac(i),i=1,ineg)
            write(*,*   ) '        Max |det J| =',jacm
          endif
          if(ineg.eq.10) then ! 10-node tet, try renumbering
            l     = ix(2)
            ix(2) = ix(1)
            ix(1) = l
            l     = ix(7)
            ix(7) = ix(6)
            ix(6) = l
            l     = ix(9)
            ix(9) = ix(8)
            ix(8) = l
          endif
        endif
      endif

2000  format(' >Element',i9,' coordinates not input for nodes:'/
     &      ('                   Local =',i3,' Global =',i9))

2001  format(' >Element',i9,' has negative jacobian at nodes:'/
     &      ('        Local =',i3,' Global =',i9,' Jacobian =',
     &       1p,1e12.4))

      end
