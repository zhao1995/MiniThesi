c$Id:$
      subroutine rrupdt(du,u,ixt,irb,numnp,ndm,ndf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Transfer rigid body incremental rotations to nodes

c      Inputs:
c         du(*)     - Increment to parameters
c         ixt(*)    - Rigid body for each node
c         irb(*)    - Rigid body equation numbers
c         numnp     - Number of nodes in mesh
c         ndm       - Spatial dimension of mesh
c         ndf       - Number of dof/node

c      Outputs:
c         u(ndf,*)  - Rotation increments in rigid body node locations
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'rigid1.h'

      integer   i,n,ndm,ndf,numnp, ixt(*),irb(nrbdof,nrbody)
      real*8    du(*), u(ndf,numnp,*)

      save

      do n = 1,numnp

c       Check for a rigid node

        if(ixt(n).gt.0) then

c         Transfer rotation increment to node

          do i = ndm+1,nrbdof
            if(irb(i,ixt(n)).gt.0) then
              u(i,n,1) = du(irb(i,ixt(n))) + u(i,n,1)
              u(i,n,2) = du(irb(i,ixt(n))) + u(i,n,2)
              u(i,n,3) = du(irb(i,ixt(n)))
            endif
          end do ! i

        end if

      end do ! n

      end
