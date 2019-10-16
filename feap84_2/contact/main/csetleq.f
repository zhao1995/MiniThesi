c$Id:$
       subroutine csetleq(ixl,nx,ilm,nl,nlag,nren,iad,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set Contact LAGrange multiplier equation numbers

c      Inputs:
c        ixl(*)   - List of nodes of contact 'element'
c        nx       - Number of nodes in 'ixl'
c        ilm(*)   - List of nodes with Lagrange multipliers
c        nl       - Number of nodes in 'ilm'
c        nlag     - Number of multipliers at each node (same at each)
c        nren(*)  - Nodal renumber array
c        numnp    - Number of nodes

c      Outputs: Passed back in global array: IAD(numnp,3)
c        IAD(i,2) - stores maximum node number for each multiplier node
c        IAD(i,3) - stores Maximum no. multipliers attached to node.
c-----[--+---------+---------+---------+---------+---------+---------+-]

       implicit   none

       integer    nx,nl,nlag, numnp, ndmax, n,nn,nr
       integer    ixl(*),ilm(*), nren(*), iad(numnp,3)

       ndmax = 0
       nr    = 0
       do n = 1,nx
         if(ixl(n).gt.0) then
           if(nren(ixl(n)).gt.ndmax) then
             ndmax = nren(ixl(n))
             nr    = ixl(n)
           endif
         endif
       end do ! n

       do n = 1,nl
         nn = ilm(n)
         if(nn.gt.0) then
           if(iad(nn,2).eq.0) then
             iad(nn,2) = nr
           elseif(ndmax.gt.nren(iad(nn,2))) then
             iad(nn,2) = nr
           endif
           iad(nn,3) = max(nlag ,iad(nn,3))
         endif
       end do ! n

       end
