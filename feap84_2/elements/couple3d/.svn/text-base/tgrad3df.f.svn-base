c$Id:$
      subroutine tgrad3df( shp3, ul, gradt,tg, ndf,nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    14/01/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute thermal gradient for 3-d elements

c      Inputs:
c        shp3(4,*)   - Shape functions and derivatives in current coord
c        ul(ndf,*)   - Nodal parameters
c        ndf         - DOF's/node (maximum)
c                      4 = temperature at nodes
c        nel         - Nodes/element

c      Outputs:
c        gradt(3)    - Thermal gradients
c        tg          - Temperature
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    ndf,nel, i,j
      real*8     shp3(4,*), ul(ndf,*), gradt(3), tg

c     Compute gradients

      do i = 1,3
        gradt(i) = 0.0d0
        do j = 1,nel
          gradt(i) = gradt(i) + shp3(i,j)*ul(4,j)
        end do ! j
      end do ! i

c     Compute temperature

      tg = 0.0d0
      do j = 1,nel
        tg = tg + shp3(4,j)*ul(4,j)
      end do ! j

      end
