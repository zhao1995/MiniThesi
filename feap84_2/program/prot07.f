c$Id:$
      subroutine prot07 (tn,ta,tl,v11,v01,v12,v02,v13,v03,dm,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c         Rotational Transformation Updates for linear problems
c         Spatial angle  3 DOF

c     INPUT variables
c         tn = 3x3 Lambda at t_n.
c         vn = 3x1 Spatial rotational velocity at t_n.
c         an = 3x1 Spatial rotational acceleration at t_n.
c         dm = 3x1 Spatial rotational increment at t_n+1 (Solver)
c         isw=     Task switch

c     OUTPUT variables
c         ta  = 3x3 Lambda at t_n+a.
c         tl  = 3x3 Lambda at t_n+1.
c         v11 = 3x1 First column of matrix velocity at time t_n+1
c         v12 = 3x1 Second column of matrix velocity at time t_n+1
c         v13 = 3x1 Third column of matrix velocity at time t_n+1
c         v01 = 3x1 First column of matrix velocity at time t_n
c         v02 = 3x1 Second column of matrix velocity at time t_n
c         v03 = 3x1 Third column of matrix velocity at time t_n

c     REMARK
c         Implemented for Enrgy-Momentum method
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'fdata.h'
      include  'tdata.h'

      integer   isw   , i   , j
      real*8    tn(4,2) , ta(4,2) , tl(4,2) , dm(3)
      real*8    v01(3)  , v02(3)  , v03(3)  , v11(3)  , v12(3) , v13(3)
      real*8    du(3)

      save

c     INITIALIZE THE ARRAYS:

      if(isw.eq.0) then

        do i = 1,2
          do j = 1,3
            ta(j,i) = 0.0d0
            tn(j,i) = 0.0d0
            tl(j,i) = 0.0d0
          end do ! j
          ta(4,i) = 1.d0
          tn(4,i) = 1.d0
          tl(4,i) = 1.d0
        end do ! i

        return
      endif

c     AT BEGINING OF A TIME STEP:

      if (isw.eq.1 .or. isw.eq.3) then

        call pzero ( du , 3 )

c     Compute spatial increment in an iteration

      else

        do i = 1 , 3
          du(i) = dm(i)
        end do ! i

      endif

      end
