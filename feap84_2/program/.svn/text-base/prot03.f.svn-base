c$Id:$
      subroutine prot04 (tn,ta,tl,du,ua,vn,v1,an,a1,dm,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Left SO(3) Rotational Transformation Updates: (Dtheta)
c               Rotational Updates using the Exponential Map

c      Inputs:
c         tn(4)    - 1-4 Quaternions at n.
c         vn(3)    - 3x1 Spatial rotational velocity at t_n.
c         an(3)    - 3x1 Spatial rotational acceleration at t_n.
c         dm(3)    - 3x1 Spatial Delta-theta at time t_n+a (solver)
c         isw      - Task switch

c      Outputs:
c         du(3)    - 3x1 Spatial rotational increment at t_n+a.
c         ta(8)    - 1-4 Quaternions at n+a.
c                    5-8 Relative Quaternions tn-tn+a.
c         tl(4)    - 1-4 Quaternions at n+1.
c         du(3)    - 3x1 Spatial Delta-theta at time t_n+a
c         ua(3)    - 3x1 Spatial rotation at t_n+a.
c         v1(3)    - 3x1 Spatial rotational velocity at t_n+1.
c         a1(3)    - 3x1 Spatial rotational acceleration at t_n+1.

c       REMARK
c         Quaternions are Stored as: (vector,scalar).
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'tdata.h'
      include  'ddata.h'
      include  'fdata.h'

      integer   isw , i, iopstr
      real*8    faca, facv , ws
      real*8    tn(9), ta(9), tl(9), dm(3), dq(4), qt(4)
      real*8    du(3), ua(3), vn(3), v1(3), an(3), a1(3)
      real*8    u1(3), q1(4), vs(3), as(3), vr(3), ar(3)

      save

      data      iopstr /0/

c     INITIALIZE QUATERNION:

      if(isw.eq.0) then

        call pzero  ( tn    , 9  )
        tn(4) = 1.d0

        call pzero  ( ta    , 9  )
        ta(4) = 1.d0
        ta(8) = 1.d0

        call pzero  ( tl    , 9  )
        tl(4) = 1.d0

        return
      endif


c     AT BEGINING OF A TIME STEP:

      if (isw.eq.1 .or. isw.eq.3) then

        if(fl(9)) then

c         Move {v,a}n+1 -> {v,a}n:

          if (isw.eq.1) then
            do i = 1,3
              vn(i) = v1(i)
              an(i) = a1(i)
            end do ! i
          else
            do i = 1,3
              v1(i) = vn(i)
              a1(i) = an(i)
            end do ! i
          endif

c         Initial Guess for Relative Velocity: tn-tn+a:

c         Zero Increment:
          if (iopstr.eq.0) then
            call pzero ( dm , 3 )

c         Zero Accelaration:
          elseif (iopstr.eq.1) then
            do i = 1 , 3
               dm(i) = theta(3)*dt*(vn(i)+(0.5d0-theta(1))*an(i)*dt)
            end do ! i

c         Increment with Previous Relative Quaternion:
          elseif (iopstr.eq.2) then
            do i = 1,3
              ua(i) = dm(i)
            end do ! i
          endif
        endif

c       Zero Relative Quaternion:

        do i = 5,7
          ta(i) = 0.0d0
        end do ! i
        ta(8) = 1.d0

      endif

c     FOR AN ITERATION WITHIN A TIME STEP (STATIC)

c     Normalize All Quaternions:

      call quanrm ( tn    )
      call quanrm ( ta    )
      call quanrm ( ta(5) )
      call quanrm ( tl    )

c     Copy dm into du and extract incremental Quaternion:

      do i = 1,3
        du(i) = dm(i)
      end do ! i
      call rotqua ( dm , dq )

c     Update relative Quaternion:

      call quamul ( dq , ta(5) , qt )
      do i = 1,4
        ta(i+4) = qt(i)
      end do ! i

c     Update total Quaternion at t-n+a and extract rotation vector

      call quamul ( ta(5) , tn(1) , qt )
      do i = 1,4
        ta(i) = qt(i)
      end do ! i
      call quarot ( ta(5) , ua )

c     Compute Relative Rotation Vector tn-tn+1

      if (theta(3) .eq. 0.d0) theta(3) = 1.d0
      do i = 1,3
        u1(i) = ua(i) / theta(3)
      end do ! i

c     Compute relative and total Quaternions at t-n+1:

      call rotqua ( u1 , q1 )
      call quamul ( q1 , tn , tl )

c     DYNAMIC UPDATES FOR AN ITERATION WITHIN A TIME STEP

      if (fl(9)) then

c       Intermediate Velocity and Acceleration Vectors:

        do i = 1 , 3
          ws    = vn(i) + 0.5d0*dt*an(i)
          vs(i) = vn(i) + (theta(2)/theta(1)-2.d0) * ws
          as(i) =(vn(i) + (0.5d0   /theta(1)-1.d0) * ws) * 2.d0/dt
        end do ! i

c       Transform Intermediate Vectors:

        call quavec ( q1 , vs , vr )
        call quavec ( q1 , as , ar )

c       Update Velocities and Accelerations:

        facv     = theta(2) / (theta(3) * theta(1) * dt)
        faca     = 1.d0     / (theta(3) * theta(1) * dt**2)
        do i = 1 , 3
          v1(i) = facv * ua(i) - vr(i)
          a1(i) = faca * ua(i) - ar(i)
        end do ! i

      endif

      end
