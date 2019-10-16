c$Id:$
      subroutine snblk(nm,ndm, ixl,xl, ns)

c      * * F E A P * * A Finite Elemen2 Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generates normal vectors for master block

c      Inputs:
c         nm        - Number of nodes on master element
c         ndm       - Spatial dimension ofmesh
c         ixl(*)    - List of nodes on master element
c         xl(ndm,*) - Nodal coordinates of master element

c      Outputs:
c         ns(3,*)   - Normals to master nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   nm,ndm, i,j,k,ii, ixl(9)
      real*8    xsj,rlen,xl(ndm,9),xs(3,2),ns(3,9),shp(3,9),ss(2,9)

      save

      data ss/-1.d0,-1.d0,  1.d0,-1.d0,  1.d0, 1.d0, -1.d0, 1.d0,
     &         0.d0,-1.d0,  1.d0, 0.d0,  0.d0, 1.d0, -1.d0, 0.d0,
     &         0.d0, 0.d0/

c     Compute normals for each bottom surface master node

      do ii = 1,nm
        if(ixl(ii).gt.0) then
          call shp2d(ss(1,ii),xl,shp,xsj,ndm,nm,ixl,.true.)

          do i = 1,ndm
            do j = 1,2
              xs(i,j) = 0.0d0
              do k = 1,nm
                xs(i,j) = xs(i,j) + xl(i,k)*shp(j,k)
              end do ! k
            end do ! j
          end do ! i

          ns(1,ii) = xs(2,1)*xs(3,2) - xs(3,1)*xs(2,2)
          ns(2,ii) = xs(3,1)*xs(1,2) - xs(1,1)*xs(3,2)
          ns(3,ii) = xs(1,1)*xs(2,2) - xs(2,1)*xs(1,2)

          rlen = 1.d0/sqrt(ns(1,ii)**2+ns(2,ii)**2+ns(3,ii)**2)

          ns(1,ii) = ns(1,ii)*rlen
          ns(2,ii) = ns(2,ii)*rlen
          ns(3,ii) = ns(3,ii)*rlen

        end if

      end do ! ii

      end
