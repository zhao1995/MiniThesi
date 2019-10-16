c$Id:$
      subroutine rcglam(rcg,rlam,nrbody)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set center of mass position and initial rotations

c      Inputs:
c         rcg(3,*,*) - Center of mass for reference state only
c         nrbody     - Number of rigid bodies

c      Outputs:
c         rcg(3,*,*) - Center of mass for initial conditions
c         rlam(9,*,*) - Initial rotation quaternion
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,n,nrbody
      real*8    rcg(3,11,*),rlam(9,6,*)

      save

      call pzero(rlam,54*nrbody)

      do n = 1, nrbody

c       Set rcg to initial values

        do i = 1,3
          rcg(i,2,n)  = rcg(i,1,n)
          rcg(i,3,n)  = rcg(i,1,n)
          rcg(i,4,n)  = 0.0d0
          rcg(i,5,n)  = 0.0d0
          rcg(i,6,n)  = 0.0d0
          rcg(i,7,n)  = 0.0d0
          rcg(i,8,n)  = 0.0d0
          rcg(i,9,n)  = 0.0d0
          rcg(i,10,n) = 0.0d0
          rcg(i,11,n) = 0.0d0
        end do ! i

c       Set rotation to inital values

        do i = 1,3
          rlam(4,i,n) = 1.d0
        end do ! i
        rlam(8,2,n) = 1.d0
      end do ! n

      end
