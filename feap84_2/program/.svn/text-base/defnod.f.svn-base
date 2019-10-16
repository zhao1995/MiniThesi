c$Id:$
      subroutine defnod(x,u,id, ld,x2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Set position of flexible node for joint and
c               assign assembly data in ld.

c     Inputs:
c        x(*)      - Reference system coordinates of node
c        id(*)     - Equation number for node
c        u(*)      - Displacements at t_n+1

c     Outputs:
c        ld(3,*)   - Joint assembly array
c        x2        - Position of node
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   j, ld(3,2), id(*)
      real*8    x(*), u(*), x2(3)

      save

      do j = 1,3
        ld(j,1) = id(j)
        ld(j,2) = 0
        x2(j)   = x(j) + u(j)
      end do ! j

      end
