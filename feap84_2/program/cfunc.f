c$Id:$
      subroutine cfunc(shp,xl,ixl,nm,ndm,x)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute coordinates for point defined by local arrays.

c      Inputs:
c         shp(3,*)  - Shape function array
c         xl(ndm,*) - Array of element coordinates
c         ixl(*)    - Element node numbers
c         nm        - Number of nodes to check
c         ndm       - Spatial dimension of mesh

c      Outputs:
c         x(ndm)    - Coordinates of point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   k,l,nm, ndm, ixl(*)
      real*8    shp(3,*),xl(3,*),x(*)

      save

      do l = 1,ndm
        x(l) = 0.0d0
      end do ! l

      do k = 1,nm
        if(ixl(k).gt.0) then
          do l = 1,ndm
            x(l) = x(l) + shp(3,ixl(k))*xl(l,ixl(k))
          end do ! l
        end if
      end do ! k

      end
