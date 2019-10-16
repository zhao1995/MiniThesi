c$Id:$
      subroutine updrot(u,ndf,xlg,mropt,numnp,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'xln' to 'xlg' to avoid conflict with     30/12/2007
c          data in 'erotas.h'
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Driver for global rotational update algorithms
c               ASSUMPTION (for now):
c                - Rotations are in partition 1
c                - Rotational DOF start at 4

c      Inputs:
c         u(ndf,numnp)  - Updated incremental director displacement
c         ndf           - Total number DOF/node
c         xlg(9,6,numnp)- Nodal rotation matrices
c         mropt(numnp,2)- Type of nodal rotational updated
c         isw           - Task (isw=1,2,3)

c      Outputs:
c         xlg(9,6,numnp)- Updated nodal rotations
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   j,isw,ndf,numnp,n,ng
      integer   mropt(numnp,2)
      real*8    u(ndf,numnp),xlg(9,6,numnp)

      save

c     Update rotation matrices to begin a step

      if (isw.eq.1) then

        do n = 1,numnp
          ng = mropt(n,2)
          if(ng.gt.0) then
            do j = 1,9
              xlg(j,1,ng) = xlg(j,3,ng)
              xlg(j,2,ng) = xlg(j,3,ng)
              xlg(j,4,ng) = xlg(j,5,ng)
            end do ! j
          endif
        end do ! n

c     Backup solution vectors to reinitiate a step

      elseif (isw.eq.3) then

        do n = 1,numnp
          ng = mropt(n,2)
          if(ng.gt.0) then
            do j = 1,9
              xlg(j,2,ng) = xlg(j,1,ng)
              xlg(j,3,ng) = xlg(j,1,ng)
              xlg(j,5,ng) = xlg(j,4,ng)
            end do ! j
          endif
        end do ! n
      endif

c     Transfer to correct update routine

      do n = 1,numnp
        ng = mropt(n,2)
        if(ng.gt.0) then

c         User rotational updates

          if (mropt(n,1).gt.0) then
            if (mropt(n,1).eq.1) then
              call urot01(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)
            elseif (mropt(n,1).eq.2) then
              call urot02(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)
            elseif (mropt(n,1).eq.3) then
              call urot03(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)
            elseif (mropt(n,1).eq.4) then
              call urot04(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)
            elseif (mropt(n,1).eq.5) then
              call urot05(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)
            elseif (mropt(n,1).eq.6) then
              call urot06(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)
            elseif (mropt(n,1).eq.7) then
              call urot07(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)
            elseif (mropt(n,1).eq.8) then
              call urot08(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)
            elseif (mropt(n,1).eq.9) then
              call urot09(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)
            elseif (mropt(n,1).eq.10) then
              call urot10(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)
            endif

c         Program rotational updates

          elseif(mropt(n,1).lt.0) then

c           Used for 'shl3df' 5-dof case

            if (mropt(n,1).eq.-1) then
              call prot01(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)

c           Used for 'framf3b' and 'framf3c'

            elseif (mropt(n,1).eq.-2) then
              call prot02(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)

c           Used for 'framf3d'

            elseif (mropt(n,1).eq.-3) then
              call prot03(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)

c           Used for 'framf3e'

            elseif (mropt(n,1).eq.-4) then
              call prot04(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)

c           Used for 'shl3df' for 6-dof

            elseif (mropt(n,1).eq.-5) then
              call prot05(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)

c           Used for explicit rigid body

            elseif (mropt(n,1).eq.-6) then
              call prot06(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)

c           Linear updates

            elseif (mropt(n,1).eq.-6) then
              call prot07(xlg(1,1,ng),xlg(1,2,ng),xlg(1,3,ng),
     &                    xlg(1,4,ng),xlg(1,5,ng),
     &                    xlg(4,4,ng),xlg(4,5,ng),
     &                    xlg(7,4,ng),xlg(7,5,ng),u(4,n),isw)
            endif
          endif
        endif
      end do ! n

      end
