c$Id:$
      subroutine pcartre(r,ang,ndf,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    22/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Convert polar components to Cartesian components for
c               outputs and plots

c      Input:
c         r(ndf,*) - Polar components of reactions
c         ang(*)   - Angle at each node
c         ndf      - Dof's/node
c         numnp    - Number of nodes

c      Outputs:
c         r(ndf,*) - Cartesian components of reactions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'corset.h'
      include   'p_point.h'

      include   'pointer.h'
      include   'comblk.h'

      integer    ndf,numnp, n
      real*8     cs,sn
      real*8     r(ndf,numnp),ang(numnp), rr(3)

c     Loop over nodes

      do n = 1,numnp
        if(anglefl) then
          if(ang(n).ne.0.0d0) then
            call pdegree(ang(n), sn,cs)
            rr(1)  = r(1,n)
            rr(2)  = r(2,n)
            r(1,n) = rr(1)*cs - rr(2)*sn
            r(2,n) = rr(1)*sn + rr(2)*cs
          endif
        elseif(eulerfl) then
        elseif(triadfl) then
          point = np(274) + 9*n - 10
          if(hr(point+1).gt.-100.0d0) then
            rr(1)  = r(1,n)
            rr(2)  = r(2,n)
            rr(3)  = r(3,n)
            r(1,n) = hr(point+1)*rr(1)
     &             + hr(point+4)*rr(2)
     &             + hr(point+7)*rr(3)
            r(2,n) = hr(point+2)*rr(1)
     &             + hr(point+5)*rr(2)
     &             + hr(point+8)*rr(3)
            r(3,n) = hr(point+3)*rr(1)
     &             + hr(point+6)*rr(2)
     &             + hr(point+9)*rr(3)
          endif
        endif
      end do ! n

      end
