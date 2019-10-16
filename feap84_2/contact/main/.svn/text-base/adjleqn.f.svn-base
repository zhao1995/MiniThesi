c$Id:$
      subroutine adjleqn(iad,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute 'offset' of each Lagrange multiplier associated
c               to 'maximum' node number (some multipliers could share
c               same maximum node number.

c      Inputs:
c        iad(n,*) - 2: stores maximum node associated with multiplier
c                      node 'n';
c                   3: stores number of multiplers for node.
c        numnp    - Number of nodes in mesh.

c      Output:
c        iad(n,*) - 1: has total number of multiplier equation of node
c                   iad(n,2).
c                 - 3: has offset for this node.
c                   N.B. multiplier equations will be placed after
c                   equations of node iad(n,2) with offset in iad(n,3).
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    numnp, n,nn
      integer    iad(numnp,*)

c     Set inital offsets to zero

      do n = 1,numnp
        iad(n,1) = 0
      end do ! n

c     Compute offset for each node.

      do n = 1,numnp
        nn = iad(n,2)                      ! Max renumbered node for 'n'
        if(nn.gt.0) then
          iad(nn,1) = iad(nn,1) + iad(n,3) ! Final equation for 'nn'
          iad(n ,3) = iad(nn,1)            ! Final equation for 'n'
        endif
      end do ! n

      end
