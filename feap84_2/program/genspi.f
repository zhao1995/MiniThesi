c$Id:$
      subroutine genspi(v,ia,ix,x,nen,nen1,ndf,ndm,numel,numnp,omg,
     &                  prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate initial velocites for spinning body

c      Inputs:
c         ix(nen1,*) - Element nodal connection list
c         x(ndm,*)   - Nodal coordinates
c         nen        - Number of nodes/element maximum
c         nen1       - Dimension for ix array
c         ndf        - Number of dof/node
c         ndm        - Spatial dimension of mesh
c         numel      - Number of elements in mesh
c         numnp      - Number of nodes in mesh
c         omg(*)     - Spin angular velocity
c         prt        - Print generated data if true
c         prth       - Print title/header data if true

c      Outputs:
c         v(ndf,*)   - Nodal velocities due to spin

c      Scratch:
c         ia(*)      - Array to store active regions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   prt,prth
      integer   i,n,ndf,nen,nen1,ndm,numel,numnp
      integer   ia(numnp),ix(nen1,numel)
      real*8    v(ndf,numnp),x(ndm,numnp),omg(3)

      save

c     Look for active regions

      do n = 1,numnp
        ia(n) = 0
      end do ! n

      do n = 1,numel
        if(ix(nen1-1,n).ge.0) then
          do i = 1,nen
            if(ix(i,n).gt.0 .and. ix(i,n).le.numnp) then
              ia(ix(i,n)) = 1
            end if
          end do ! i
        endif
      end do ! n

c     Set velocities at nodes

      if(ndm.eq.2 .and. ndf.ge.2) then
        do n = 1,numnp
          if(mr(np(190)+n-1).ge. 0 .and. ia(n).gt.0) then
            v(1,n) =  x(2,n)*omg(3)
            v(2,n) = -x(1,n)*omg(3)
          endif
        end do ! n
      elseif(ndm.eq.3 .and. ndf.ge.3) then
        do n = 1,numnp
          if(mr(np(190)+n-1).ge. 0 .and. ia(n).gt.0) then
            v(1,n) =  x(2,n)*omg(3) - x(3,n)*omg(2)
            v(2,n) =  x(3,n)*omg(1) - x(1,n)*omg(3)
            v(3,n) =  x(1,n)*omg(2) - x(2,n)*omg(1)
          endif
        end do ! n

      endif

c     Output velocities

      if(prt) then
        call prtitl(prth)
        write(iow,2000) (i,omg(i),i=1,3),(i,i=1,ndm)
        if(ior.lt.0) then
          write(*,2000) (i,omg(i),i=1,3),(i,i=1,ndm)
        endif
        do n = 1,numnp
          if(mr(np(190)+n-1) .ge. 0) then
            write(iow,2001) n,(v(i,n),i=1,ndm)
            if(ior.lt.0) then
              write(*,2001) n,(v(i,n),i=1,ndm)
            endif
          endif
        enddo ! n
      endif

2000  format('   I n i t i a l   V e l o c i t i e s'//'  Spin:',
     &     3(' omg(',i1,') = ',1p,e11.3)/'    Node',3(i6,'-velo'))

2001  format(i8,1p,3e11.3)

      end
