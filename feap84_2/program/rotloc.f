c$Id:$
      subroutine rotloc(xlg,tg,mropt,numnp,nl,ngn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'xg' to 'xlg' (same as in updrot          30/12/2007
c       2. change 'ng' to 'ndn' in set of ndof from mropt   01/11/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Localize rotational arrays for shells.

c      Inputs:
c         xlg(9,6,*)  - Global rotation matrix
c         tg(*)       - Global fiber thickness
c         mropt(*,2)  - Nodal rotation update option
c         numnp       - Number of nodes in mesh
c         nl          - Local node number
c         ngn         - Global node number

c      Outputs:
c         xln(9,*,*)  - Local rotation matrix
c         rots(3,*,*) - Local rotations
c         rvel(3,*,*) - Local rotation velocity
c         racc(3,*,*) - Local rotation acceleration
c         thkl(*)     - Local fiber thickness
c         ndof        - Local dof map for 5/6 rotations
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'crotas.h'
      include  'erotas.h'

      integer   numnp, nl,ngn, ng,nn
      integer   mropt(numnp,*),ndofop(10)
      real*8    xlg(9,6,numnp),tg(numnp)

      save

      data      ndofop /5,6,6,6,6,6,6,6,6,6/

      ng = mropt(ngn,2)

c     Localize rotation matrix

      if(ng.gt.0) then
        call pmove(xlg(1,1,ng),xln(1,nl,1),9)
        call pmove(xlg(1,2,ng),xln(1,nl,2),9)
        call pmove(xlg(1,3,ng),xln(1,nl,3),9)
        call pmove(xlg(1,6,ng),xln(1,nl,4),9)

c       Localize rotations

        call pmove(xlg(1,4,ng),rots(1,nl,1),3)
        call pmove(xlg(1,5,ng),rots(1,nl,2),3)

c       Localize rotational velocities

        call pmove(xlg(4,4,ng),rvel(1,nl,1),3)
        call pmove(xlg(4,5,ng),rvel(1,nl,2),3)

c       Localize rotational accelerations

        call pmove(xlg(7,4,ng),racc(1,nl,1),3)
        call pmove(xlg(7,5,ng),racc(1,nl,2),3)

c       Localize fiber thickness

        thkl(nl) = tg(ng)

c       Localize dof map for 5 and 6 d.o.f.'s

        if(nl.le.9) then
          nn   = -mropt(ngn,1)
          if(nn.gt.0 .and. nn.lt.10) then
            ndof(nl) = ndofop(nn)
          else
            ndof(nl) = 0
          endif
        endif
      endif

      end
