c$Id:$
      subroutine findxy(ixy,sne,x,ix,nen,scale,s0,sx, flge)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Finds intersection points with elements

c      Inputs:
c         ixy(2,*)  - Screen coordinate array
c         x(3,*)    - Nodal coordinates of element
c         ix(*)     - Nodal numbers for element
c         scale     - Scale of plot
c         s0(2)     - Center of plot coordinates
c         sx(*)     - Origin of plot


c      Outputs:
c         sne(2,*)  - Intersection data for element
c         flge      - Flag, element has intersections if true
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata2.h'

      logical   flge,plus,negn
      integer   nen,k,n
      real*8    scale,l,l2,dn,x1,x2,y1,y2,ddx,ddy,s1,s2

      integer   ix(nen),ixy(2,2)
      real*8    x(3,*),sne(2,nen),s0(2),sx(2)

      save

c     Compute coordinates and differences

      x1  = ixy(1,1)/dble(idx)
      x2  = ixy(1,2)/dble(idx)
      y1  = ixy(2,1)/dble(idy)
      y2  = ixy(2,2)/dble(idy)
      ddx  = x2 - x1
      ddy  = y2 - y1
      l2   = ddx*ddx + ddy*ddy
      plus = .false.
      negn = .false.

      do k = 1,nen
        n = abs(ix(k))
        if(n.gt.0) then

c         Compute normal coordinates

          s1 = scale*(x(1,n) + x(1,n) - sx(1)) + s0(1)
          s2 = scale*(x(2,n) + x(2,n) - sx(2)) + s0(2)

c         Check if line passes through element

          l = (s1-x1)*ddx + (s2-y1)*ddy
          if(l.gt.0.0d0 .and. l.lt.l2) then
            dn = (x1-s1)*ddy + (s2-y1)*ddx
            if(dn.ge.0.0d0) plus = .true.
            if(dn.le.0.0d0) negn = .true.
          endif
          if(flge) then
            sne(1,k) = l/l2
            sne(2,k) = dn/l2
          endif
        endif
      end do ! k

      flge = plus .and. negn

      end
