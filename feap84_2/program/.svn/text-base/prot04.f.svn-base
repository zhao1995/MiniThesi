c$Id:$
      subroutine prot03 (tn,tl,t1,du,u1,vn,v1,an,a1,dm,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Left SO(3) Rotational Transformation Updates: (Dtheta)
c               Rotational Updates using Cayley Transformation
c               Used with energy-momentum updates (noi=5).

c      Inputs:
c         tn(4)    - 1-4 Quaternion at n.
c         vn(3)    - 3x1 Spatial rotational velocity at t_n.
c         an(3)    - 3x1 Spatial rotational acceleration at t_n.
c         dm(3)    - 3x1 Spatial Delta-theta at time t_n+1 (solver)
c         isw      - Task switch

c      Outputs:
c         tl(8)    - 1-4 Quaternion at n+1.
c                    5-8 Relative Quaternion tn-tn+1.
c         t1(4)    - 1-4 Quaternion at n+1.
c         du(3)    - 3x1 Spatial rotational increment at time t_n+1
c         u1(3)    - 3x1 Spatial rotation at t_n+1.
c         v1(3)    - 3x1 Spatial rotational velocity at t_n+1.
c         a1(3)    - 3x1 Spatial rotational acceleration at t_n+1.

c       REMARK
c         Quaternion are Stored as: (vector,scalar).
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'tdata.h'
      include  'ddata.h'
      include  'fdata.h'

      integer   isw , i
      real*8    faca, facv , ws
      real*8    tn(9), tl(9), dm(3), dq(4), qt(4), t1(9)
      real*8    du(3), vn(3), v1(3), an(3), a1(3)
      real*8    u1(3), vs(3), as(3), vr(3), ar(3)

      save

c     INITIALIZE QUATERNION:

      if(isw.eq.0) then
        call pzero  ( tn    , 9  )
        tn(4) = 1.d0
        tn(8) = 1.d0

        call pzero  ( tl    , 9  )
        tl(4) = 1.d0
        tl(8) = 1.d0

        call pzero  ( t1    , 9  )
        t1(4) = 1.d0
        t1(8) = 1.d0

        return
      endif

c     AT BEGINING OF TIME STEP:

      if (fl(9) .and. (isw.eq.1 .or. isw.eq.3)) then

c       Move {v,a}n+1 -> {v,a}n:

        if (isw.eq.1) then
          call pmove ( v1 , vn , 3 )
          call pmove ( a1 , an , 3 )
        else
          call pmove ( vn , v1 , 3 )
          call pmove ( an , a1 , 3 )
        endif

c       Normalize initial Quaternion:

        call quanrm ( tn )

c       Initial Guess for Relative Velocity: tn-tn+1:

c       Zero Increment:

        call pzero ( dm , 3 )

      endif

      if(isw.eq.1)then

c     Zero Relative Quaternion:

        call pzero ( tl(5) , 3 )
        tl(8) = 1.d0

      endif

c     FOR ITERATION WITHIN TIME STEP (STATIC)

c     Normalize All Quaternion:

      call quanrm ( tn    )
      call quanrm ( tl    )
      call quanrm ( tl(5) )

c     Copy dm into du and extract incremental Quaternion:

      call pmove ( dm , du , 3 )
      call rqcay ( dm , dq )

c     Update relative Quaternion:

      call quamul ( dq , tl(5) , qt )
      call pmove  ( qt , tl(5) , 4  )

c     Update total Quaternion at t-n+l and extract rotation vector

      call quamul ( tl(5) , tn(1) , qt )
      call pmove  ( qt    , tl(1) , 4  )
      call qrcay  ( tl(5) , u1 )
      call pmove  ( tl    , t1    , 4  )

c     DYNAMIC UPDATES FOR ITERATION WITHIN TIME STEP

      if (fl(9)) then

c       Intermediate Velocity and Acceleration Vectors:

        do i = 1 , 3
          ws    = vn(i) + 0.5d0*dt*an(i)
          vs(i) = vn(i) + (theta(2)/theta(1)-2.d0) * ws
          as(i) =(vn(i) + (0.5d0   /theta(1)-1.d0) * ws) * 2.d0/dt
        end do ! i

c       Transform Intermediate Vectors:

        call quavec ( tl(5) , vs , vr )
        call quavec ( tl(5) , as , ar )

c       Update Velocities and Accelerations:

        facv     = theta(2) / (theta(1) * dt)
        faca     = 1.d0     / (theta(1) * dt**2)
        do i = 1 , 3
          v1(i) = facv * u1(i) - vr(i)
          a1(i) = faca * u1(i) - ar(i)
        end do ! i

      endif

      end
