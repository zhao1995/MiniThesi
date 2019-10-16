c$Id:$
      subroutine pnorm2d(enorm,xl,ndm,nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute normal vectors (unscaled) for 2-d triangular
c               and quadrilateral elements

c      Inputs:
c         xl(ndm,*)  - Element nodal coordinates
c         ndm        - Spatial dimension of mesh
c         nel        - Number of nodes on element

c      Outputs:
c         enorm(3,*) - Unscaled normal vector components
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    ndm,nel, i,j,k,nm
      real*8     enorm(3,*), xl(ndm,*)

      save

      if(nel.eq.6) then
        nm = 3
      else
        nm = min(4,nel)
      endif
      do i = 1,nm
        j = mod(i     ,nm) + 1
        k = mod(i+nm-2,nm) + 1
        enorm(1,i) = xl(2,j) - xl(2,k)
        enorm(2,i) = xl(1,k) - xl(1,j)
      end do ! i

c     Triangle: 6 node

      if(nel.eq.6) then
        do i = 1,3
          enorm(1,i+3) = 1.d0
        end do ! i

c     Quadrilateral: 8 or 9 node

      elseif(nel.eq.8 .or. nel.eq.9) then
        do i = 1,4
          enorm(1,i+4) = 1.d0
        end do ! i
        if(nel.eq.9) enorm(1,9) = 1.d0

c     Quadrilateral: 16 node

      elseif(nel.eq.16) then
        do i = 1,4
          j = mod(i,4) + 1
          enorm(1,2*i+3) = 1.d0
          enorm(1,2*i+4) = 1.d0
        end do ! i
      endif

      end
