c$Id:$
      subroutine rxnod(xh,x1,rh,r1,yh,y1,x,rcg,rlam,alpr,ixt,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Computation of rigid body position vectors

c      Inputs:
c         x(ndm,*)  - Nodal coordinates for mesh
c         rcg(3,*)  - Center of mass solution state
c         rlam(9,*) - Rotational state
c         alpr      - Parameter to obtain solution at t_n
c         ixt(*)    - Nodal rigid body numbers
c         ndm       - Spatial dimension of mesh

c      Outputs:
c         xh(3)     - Node position at mid-step
c         x1(3)     - Node position at t_n+1
c         rh(3)     - Center of mass position at mid-step
c         r1(3)     - Center of mass position at t_n+1
c         yh(3)     - Center of mass to node postion at t_n+1/2
c         y1(3)     - Center of mass to node postion at t_n+1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, ixt, ndm
      real*8    alpr,dx(3),x(3),  rn(3),rh(3),r1(3)
      real*8    xn(3),xh(3),x1(3),yn(3),yh(3),y1(3)
      real*8    rcg(3,11,*),rlam(9,6,*)

      save

c     Compute positions at t_n and t_n+1

      do i = 1,ndm
        dx(i) = x(i) - rcg(i,1,ixt)
      end do ! i

      call quavec(rlam(1,1,ixt),dx,yn)
      call quavec(rlam(1,3,ixt),dx,y1)

      do i = 1,ndm
        r1(i) = rcg(i,2,ixt)
        rn(i) = r1(i) - rcg(i,3,ixt)*alpr
        x1(i) = r1(i) + y1(i)
        xn(i) = rn(i) + yn(i)
      end do ! i

      do i = 1,ndm
        xh(i) = (xn(i) + x1(i))*0.5d0
        rh(i) = (rn(i) + r1(i))*0.5d0
        yh(i) = (yn(i) + y1(i))*0.5d0
      end do ! i

      end
