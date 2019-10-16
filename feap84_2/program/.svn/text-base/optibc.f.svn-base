c$Id:$
      subroutine optibc(id,nnid,ndf,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Count active dof at each node

c     Inputs:
c       ndf        - Number dof at each node
c       numnp      - Number nodes
c       id(ndf,*)  - Equation numbers for nodes

c     Outputs:
c       nnid(*)    - Number active dof/node
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    ndf,numnp, i,n
      integer    id(ndf,numnp),nnid(numnp)

      save

      do n = 1,numnp
        nnid(n) = 1
        nnid(n) = 0
        do i = 1,ndf
          if(id(i,n).gt.0) then
            nnid(n) = nnid(n) + 1
          endif
        end do ! i
      end do ! n

      end
