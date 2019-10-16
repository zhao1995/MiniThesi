c$Id:$
      subroutine prot05 (tn,ta,tl,v11,v01,v12,v02,v13,v03,dm,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c         Rotational Transformation Updates for Shell Intersection
c         Spatial Delta-theta  3 DOF
c         Exponential updates in terms of matrices, no quaternions
c         DIRECTOR VELOCITIES

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

c     IMPORTANT variable routines
c         dl = 3x3 Delta-Lambda.
c     REMARK
c         Implemented for Enrgy-Momentum method
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'fdata.h'
      include  'tdata.h'

      integer   isw   , i   , j
      real*8    dtr     , tn(3,3) , ta(3,3) , tl(3,3) , dl(3,3) , dm(3)
      real*8    v01(3)  , v02(3)  , v03(3)  , v11(3)  , v12(3) , v13(3)
      real*8    tt(3,3) , du(3)

      save

c     INITIALIZE THE ARRAYS:

      if(isw.eq.0) then

        do i = 1,3
          do j = 1,3
            ta(j,i) = 0.0d0
            tn(j,i) = 0.0d0
            tl(j,i) = 0.0d0
          end do ! j
          ta(i,i) = 1.d0
          tn(i,i) = 1.d0
          tl(i,i) = 1.d0
        end do ! i

        return
      endif

c     AT BEGINING OF A TIME STEP:

      if (fl(9) .and. (isw.eq.1 .or. isw.eq.3)) then

c     Move {v,a}n+1 -> {v,a}n

        if (isw.eq.1) then
          call pmove ( v11 , v01 , 3 )
          call pmove ( v12 , v02 , 3 )
          call pmove ( v13 , v03 , 3 )
        else
          call pmove ( v01 , v11 , 3 )
          call pmove ( v02 , v12 , 3 )
          call pmove ( v03 , v13 , 3 )
        endif

        call pzero ( du , 3 )

c     Compute spatial increment in an iteration

      else

        do i = 1 , 3
          du(i) = dm(i)
        end do ! i

      endif

c     UPDATE WITHIN TIME STEP

c     Compute rotation matrix and quaternion for du

      call lamrot(du,dl)

c     Update rotation matrix at t_n+a

      call pmove ( tl , tt , 9 )

      do i = 1 , 3
        do j = 1 , 3
           tl(i,j) = dl(i,1)*tt(1,j) + dl(i,2)*tt(2,j)
     &             + dl(i,3)*tt(3,j)
        end do ! j
      end do ! i

c     Update rotation matrix at t_n+alpha

      if(fl(9)) then
        call pmove ( ta , tt , 9 )
        do i=1,3
          du(i) = du(i)*theta(3)
        end do ! i
        call lamrot(du,dl)
        do i = 1 , 3
          do j = 1 , 3
            ta(i,j) = dl(i,1)*tt(1,j) + dl(i,2)*tt(2,j)
     &              + dl(i,3)*tt(3,j)
          end do ! j
        end do ! i
      else
        call pmove (tl, ta, 9)
      endif

c     DYNAMICS: Update Lambda Velocities:

      if (fl(9)) then
        dtr = 2.d0/dt
        do i = 1,3
          v11(i)= dtr*(tl(i,1) - tn(i,1)) - v01(i)
          v12(i)= dtr*(tl(i,2) - tn(i,2)) - v02(i)
          v13(i)= dtr*(tl(i,3) - tn(i,3)) - v03(i)
        end do ! i
      endif

      end
