c$Id:$
      subroutine preaout(id,r, rl, ur, ndf,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/12/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output reaction forces to file

c      Inputs:
c         id(ndf,numnp,2) - Boundary conditions
c         r(ndf,numnp)    - Reaction values
c         ur              - Tolerance value for output filtering
c         ndf             - DOF's at node
c         numnp           - Number nodes

c      Temporary array
c         rl(ndf)         - Store printed values for nodes

c      Outputs:
c         none            - Reactions to file 'ios'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iodata.h'

      logical    idfl
      integer    ndf, numnp, i, j
      real*8     etab, ur
      integer    id(ndf,numnp,2)
      real*8     r(ndf,numnp), rl(ndf)

c     Output reaction forces to file

      do i = 1,numnp
        etab = 0.0d0
        do j = 1,ndf
          etab = max(etab,abs(r(j,i)))
        end do ! j
        if(etab.gt.ur) then              ! Filter 'zero' values
          idfl = .false.
          do j = 1,ndf
            rl(j) = -r(j,i)              ! Reverse sign to get force
            if(abs(rl(j)).lt.ur) rl(j) = 0.0d0
            if(id(j,i,2).ne.0) then      ! Set flag to mark reaction
              idfl = .true.
            else                         ! Force degree of freedom
              rl(j) = 0.0d0
            endif
          end do ! j
          if(idfl) write(ios,*) i,(rl(j),j=1,ndf)  ! Do write
        endif
      end do ! i

      end
