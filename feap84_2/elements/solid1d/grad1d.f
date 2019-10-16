c$Id:$
      subroutine grad1d(ul,shp,epsv,ndf,ndm,nel,gradu,type)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension shp(2,*)                               02/04/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: J-integral and Matrial Force Calculations for linear
c              elements

c      Inputs:
c         ul(ndf,*)  - Current nodal solution parameters
c         shp(2,*)   - Shape functions
c         epsv(*)    - Volumetric strain measure
c         ndf        - Degree of freedoms/node
c         ndm        - Mesh coordinate dimension
c         nel        - Nodes on element

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    ndf,ndm,nel, i,j,k, type
      real*8     rr,epsv
      real*8     ul(ndf,*), shp(2,*), gradu(3,3)

c     Compute displacement gradient

      do i = 1,3
        do j = 1,3
          gradu(j,i) = 0.0d0
        end do ! j
      end do ! i
      do i = 1,ndm
        do j = 1,ndm
          do k = 1,nel
            gradu(j,i) = gradu(j,i) + ul(j,k)*shp(i,k)
          end do ! k
        end do ! j
      end do ! i

      if(type.eq.3) then
        rr         = 0.0d0
        gradu(3,3) = 0.0d0
        do k = 1,nel
          rr         = rr         + ul(1,k)*shp(2,k)
          gradu(3,3) = gradu(3,3) + ul(1,k)*shp(2,k)
        end do ! k

        if(rr.gt.0.0d0) then
          gradu(3,3) = gradu(3,3)/rr
        else
          gradu(3,3) = gradu(1,1)
        endif
      endif

c     Correct gradient for mixed model

      rr = (epsv - gradu(1,1) - gradu(2,2) - gradu(3,3))/3.d0

      gradu(1,1) = gradu(1,1) + rr
      gradu(2,2) = gradu(2,2) + rr
      gradu(3,3) = gradu(3,3) + rr

      end
