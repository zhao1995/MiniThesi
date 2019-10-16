c$Id:$
      subroutine udynam(du,u,ud,nneq,ndf,ndfp,ndfo,npart,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform updates to solution vectors

c      Inputs :
c        du(*)      - Increments to solution (isw=2 only)
c        u(nneq,2)  - Solution states at t_n+1
c        ud(nneq,*) - User controlled vectors
c        nneq       - Number nodal parameters (ndf*numnp)
c        ndf        - Number degree of freedoms/node
c        ndfp(*)    - Partition number for degree of freedoms
c        ndfo(*)    - Order for number for degree of freedoms
c        npart      - Active partition number
c        isw        - Switch: 1 - initial updates at start of step
c                             2 - iterative updates withini step
c                             3 - restore solution to values at start
c                                 of step

c      Outputs:
c        ud(nneq,*) - User controlled vectors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,n,ndf,nneq,npart,isw
      integer   ndfp(*),ndfo(*)
      real*8    du(*),u(nneq,2),ud(nneq,*)

c     Loop structure for partitions

      do i = 1,ndf
        if(ndfp(i).eq.npart) then
          do n = i,nneq,ndf
c           Perform steps here
          end do
        endif
      end do ! i

      end
