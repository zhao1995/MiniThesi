c$Id:$
      subroutine pdefm(x,b,c,ndm,ndf,numnp, dr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute deformed position of nodes

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         b(ndf,*)  - Solution vector to add to coordinates
c         c         - Scale factor for added solution
c         ndm       - Dimension of x array
c         ndf       - Number dof/node
c         numnp     - Number of nodes in mesh

c      Outputs:
c         dr(3,*)   - Deformed coordinates
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdatay.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   ndm,ndf,numnp, i,n
      real*8    x(ndm,*),b(ndf,*),uu(3), dr(3,*), c

      save

      do n = 1,numnp
        if(mr(npty+n-1).ge.0) then
          do i = 1,3
            if(pdf(i).gt.0 .and. pdf(i).le.ndf) then
              uu(i) = c*b(pdf(i),n)
            else
              uu(i) = 0.0d0
            endif
          end do ! i
          if(torsfl .and. ndf.ge.3 .and. ndm.eq.2) then
            dr(1,n) = (x(1,n) + uu(1))*cos(uu(3))
            dr(2,n) = (x(2,n) + uu(2))
            dr(3,n) = (x(1,n) + uu(1))*sin(uu(3))
          else
            do i = 1,ndm
              dr(i,n) = x(i,n) + uu(i)
            end do ! i
            do i = ndm+1,3
              dr(i,n) = uu(i)
            end do ! i
          endif
        endif
      end do ! n

      end
