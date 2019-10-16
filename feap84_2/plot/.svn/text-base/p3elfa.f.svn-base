c$Id:$
      subroutine p3elfa(x,elface,norm,iface,nface,ndm,ifc1,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine number of faces attached to each node and
c               compute the normal to each face

c      Inputs:
c         x(ndm,*)      - Nodal coordinates of mesh
c         iface(ifc1,*) - Face nodes and material set data
c         nface         - Number of faces
c         ndm           - Dimension of x array
c         numnp         - Number of nodes in mesh

c      Outputs:
c         elface(*)     - Number of faces connected to each node
c         norm(3,*)     - Normal at each face
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm,nface,ifc1,numnp,ii,jj,i,n
      integer   iface(ifc1,nface),elface(*)
      real*8    x(ndm,numnp),norm(3,nface),dx1(3),dx2(3)

      save

      do n = 1,numnp + 1
        elface(n) = 0
      end do ! n

      do n = 1,nface
        if(iface(ifc1,n).gt.0) then
          do i = 1,4
            if(iface(i,n).gt.0) then
              elface(iface(i,n)) = elface(iface(i,n)) + 1
            end if
          end do ! i
          do i = 1,3
            dx1(i) = x(i,iface(3,n)) - x(i,iface(1,n))
          end do ! i
          if(iface(4,n).gt.0) then
            do i = 1,3
              dx2(i) = x(i,iface(4,n)) - x(i,iface(2,n))
            end do ! i
          else
            do i = 1,3
              dx2(i) = x(i,iface(3,n)) - x(i,iface(2,n))
            end do ! i
          endif
          norm(1,n) = dx1(2)*dx2(3) - dx1(3)*dx2(2)
          norm(2,n) = dx1(3)*dx2(1) - dx1(1)*dx2(3)
          norm(3,n) = dx1(1)*dx2(2) - dx1(2)*dx2(1)
          dx1(1)  = sqrt(norm(1,n)*norm(1,n) + norm(2,n)*norm(2,n)
     &                 + norm(3,n)*norm(3,n))
          if(dx1(1).gt.0.0d0) then
            dx1(1) = 1.0d0/dx1(1)
          else
            dx1(1) = 1.0d0
          endif
          do i = 1,3
            norm(i,n) = norm(i,n)*dx1(1)
          end do ! i
        endif
      end do ! n

      jj        = elface(1)
      elface(1) = 0
      do n = 1,numnp
        ii          = elface(n+1)
        elface(n+1) = elface(n) + jj
        jj          = ii
      end do ! n

      end
