c$Id:$
      subroutine rmascg(rmas,rcg,rinr,xcg,ma,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Finish computation for center of mass location

c      Inputs:
c         rmas(*)     - Mass of materials
c         rcg(3,11,*) - Integral of density x position
c         rinr(9,*)   - Inertia with respect to referenc point
c         xcg(*)      - Reference point for calculation
c         ma          - Material number
c         ndm         - Spatial dimension of mesh

c      Outputs:
c         rcg(3,11,*) - Location of center of mass
c         rinr(9,*)   - Inertia at center of mass
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ma, ndm, i,j
      real*8    rmcg,drcgl
      real*8    rmas(*),rinr(3,3,*),rcg(3,11,*),rcgl(3),xcg(3)

      save

      do i = 1,3
        rcgl(i) = 0.0d0
      end do ! i

      do i = 1,ndm
        if(rmas(ma).ne.0.0d0) then
          rcgl(i)     = rcg(i,1,ma)/rmas(ma)
        else
          rcgl(i) = 0.0d0
        endif
        rcg(i,1,ma) = rcgl(i) + xcg(i)
        rcg(i,2,ma) = rcg(i,1,ma)
      end do ! i

c     Adjust inertia tensor to center of mass

      if(ndm.eq.2) then
        rinr(3,3,ma) = rinr(3,3,ma)
     &               - rmas(ma)*(rcgl(1)*rcgl(1) + rcgl(2)*rcgl(2))
      elseif(ndm.eq.3) then
        drcgl = rmas(ma)*(rcgl(1)*rcgl(1) + rcgl(2)*rcgl(2)
     &                                    + rcgl(3)*rcgl(3))
        do i = 1,ndm
          rmcg = rmas(ma)*rcgl(i)
          rinr(i,i,ma) = rinr(i,i,ma) - drcgl
          do j = 1,ndm
            rinr(i,j,ma) = rinr(i,j,ma) + rmcg*rcgl(j)
          end do ! j
        end do ! i
      end if

      end
