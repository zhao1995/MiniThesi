c$Id:$
      subroutine rupjnt(du,jnt,ujnt,cc1,numjts,ndf,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Update joint solution values

c      Inputs:
c         jnt(6,*) - Joint identifiers
c         cc1      - Solution update parameter
c         numjts   - Number of joints
c         ndf      - Number dof/node
c         isw      - Switch: = 1: Initialize for new time step
c                            = 2: Increment within a time step
c                            = 3: Back-up time step

c      Outputs:
c         ujnt(3,*)- Joint parameter solutions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pointer.h'
      include  'comblk.h'

      integer   i,isw, j, n,numjts,ndf,nj
      integer   jnt(6,numjts)
      real*8    cc1, du(*),ujnt(ndf,numjts,3)

      save

c     Set initial values of increments to zero

      if(isw.eq.1) then

        do n = 1,numjts
          do i = 1,ndf
c           ujnt(i,n,1) = ujnt(i,n,1) + ujnt(i,n,2)/cc1
            ujnt(i,n,2) = 0.0d0
            ujnt(i,n,3) = 0.0d0
          end do
        end do

c     Set update within step

      elseif(isw.eq.2) then

        do n = 1,numjts
          nj = jnt(5,n)
          j  = jnt(6,n)
          if(j.gt.0) then
            do i = 1,nj
              ujnt(i,n,1) = ujnt(i,n,1) + cc1*du(j+i)
              ujnt(i,n,2) = ujnt(i,n,2) + cc1*du(j+i)
              ujnt(i,n,3) =               cc1*du(j+i)
            end do
          end if
        end do

c     Back up step

      elseif(isw.eq.3) then

        do n = 1,numjts
          do i = 1,ndf
            ujnt(i,n,1) = ujnt(i,n,1) - ujnt(i,n,2)
            ujnt(i,n,2) = 0.0d0
            ujnt(i,n,3) = 0.0d0
          end do
        end do

      end if

c     Check if explicit updates

c     call  xrujnt(jnt,mr(np(96)),hr(np(95)),hr(np(104)),
c    &             hr(np(102)),isw)

      end
