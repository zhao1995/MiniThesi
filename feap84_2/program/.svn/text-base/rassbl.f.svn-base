c$Id:$
      subroutine rassbl(rmas,rcg,rinr,emas,ecg,einr,ma,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Assemble rigid mass arrays

c      Inputs:
c         emas        - Element mass
c         ecg(3)      - Element center of mass
c         e(3,3)      - Element intertia tensor
c         ma          - Rigid body number
c         ndm         - Spatial dimension of mesh

c      Outputs:
c         rmas(*)     - Mass of rigid bodies
c         rcg(3,11,*) - Location of center of mass in reference state
c         rinr(3,3,*) - Inertia tensor in reference state
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ma, ndm, i,j
      real*8    rmas(ma),rcg(3,11,ma),rinr(3,3,ma),emas,ecg(3),einr(3,3)

      save

      rmas(ma) = rmas(ma) + emas
      do i = 1,ndm
        rcg(i,1,ma) = rcg(i,1,ma) + ecg(i)
      end do ! i
      if(ndm.eq.2) then
        rinr(3,3,ma) = rinr(3,3,ma) + einr(3,3)
      elseif(ndm.eq.3) then
        do i = 1,3
          do j = 1,3
            rinr(i,j,ma) = rinr(i,j,ma) + einr(i,j)
          end do ! j
        end do ! i
      end if

      end
