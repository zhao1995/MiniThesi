c$Id:$
      subroutine prot02 (tn, tl, t1, du, u1, vn, v1, an, a1, dm, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove 'updateopt'                               16/04/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Left SO(3) Rotational Transformation Updates:
c               Rotational Updates using exponential map
c               Dynamic updates for genl alpha method (HHT & Newark)
c               ... following Simo & Vu- Quoc 1988

c      I. Romero April 2002 (modified original Simo/Taylor urot03.f)

c      Inputs:
c         tn(4)    - 1-4 Quaternion at n.
c         vn(3)    - 3x1 spatial rotational velocity at t_n.
c         an(3)    - 3x1 spatial rotational acceleration at t_n.
c         dm(3)    - 3x1 Spatial Delta-theta at time t_n+1 (solver)
c         isw      - Task switch

c      Outputs:
c         tl(8)    - 1-4 Quaternion at n+1.
c                    5-8 Relative Quaternion tn-tn+1.
c         t1(4)    - 1-4 Quaternion at n+1.
c         du(3)    - 3x1 Spatial rotational increment at time t_n+1
c         u1(3)    - 3x1 spatial rotation at t_n+1.
c         v1(3)    - 3x1 spatial rotational velocity at t_n+1.
c         a1(3)    - 3x1 Spatial rotational acceleration at t_n+1.

c       REMARKS
c         Quaternion are Stored as: (vector,scalar).

c         The solver variable dm is incremental rotation vector, in such
c         a way that rotation Lambda, at time t_n+1 in iteration k+1 is:

c                    Lambda_n+1^(k+1) = exp[dm] Lambda_n+1^(k) , or,
c                                     = exp[u1] Lambda_n

c         Tested for feap version 7.4d
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      implicit   none

      include   'tdata.h'
      include   'ddata.h'
      include   'fdata.h'

      integer    isw , i
      real*8     faca, facv , ws

      real*8     tn(9), tl(9), dm(3), dq(4), qt(4), t1(9)
      real*8     du(3), vn(3), v1(3), an(3), a1(3)
      real*8     u1(3), vs(3), as(3), vr(3), ar(3)
      real*8     beta , gamma

      save

c     Parameters for generalized alpha method

      beta  = theta(1)
      gamma = theta(2)

c     INITIALIZE ROTATIONS

      if(isw.eq.0) then

        call pzero  ( tn    , 9  )
        tn(4) = 1.d0

        call pzero  ( tl    , 9  )
        tl(4) = 1.d0
        tl(8) = 1.d0

        call pzero  ( t1    , 9  )
        t1(4) = 1.d0

        return
      endif

c     AT BEGINING OF TIME STEP OR BACK UP A STEP

      if (isw.eq.1 .or. isw.eq.3) then

c       Dynamic updates

        if (fl(9)) then

c         Move {v,a}n+1 -> {v,a}n:

          if (isw.eq.1) then
            call pmove ( v1 , vn , 3 )
            call pmove ( a1 , an , 3 )
          else
            call pmove ( vn , v1 , 3 )
            call pmove ( an , a1 , 3 )
          endif

c         Normalize initial Quaternion:

          call quanrm ( tn )

c         Initial Guess for Relative Velocity: tn-tn+1:

c         Predictor with zero Incremental rotation.

          do i = 1,3
            dm(i) = 0.0d0
          end do ! i

        end if ! end of dynamic updates

c       Set relative rotation to zero:

        do i = 5,7
          tl(i) = 0.0d0
        end do ! i
        tl(8) = 1.d0

c       Move rotation_n+1 to rotation_n or viceversa

        if (isw.eq.1) then
          call pmove ( tl(1) , tn(1) , 4 )
        else
          call pmove ( tn(1) , tl(1) , 4 )
        endif

      endif

c     FOR ITERATION WITHIN TIME STEP (STATIC)

c     Normalize rotation at tn, at tn+1 and relative rotation tn->tn+1

      call quanrm ( tn    )
      call quanrm ( tl    )
      call quanrm ( tl(5) )

c     Copy dm (solver variable) into du (Delta theta) and extract
c     incremental Quaternion, that is, rotation (in quaternion format)
c     exp[dm]

      call pmove ( dm , du , 3 )
      call rotqua( dm , dq)

c     Update relative rotation (from Lambda_n to Lambda_n+1) using
c     exponential map.
c     ...  exp(Delta theta) Delta Lambda(k) -> Delta Lambda(k+1)

      call quamul ( dq      , tl(5)          , qt )
      call pmove  ( qt      , tl(5)          , 4  )

c     Update total Quaternion at t-n+l
c     ...      Delta Lambda(k+1)  Lambda_n = Lambda_n+1

      call quamul ( tl(5)        , tn(1) , qt )
      call pmove  ( qt           , tl(1) , 4  )

c     Extract spatial rotation vector u1 = exp^-1( Delta Lambda(k+1) ).
c     Thus, Lambda_n+1^(k+1) = exp[u1] Lambda_n
c                            = exp[dm] Lambda_n+1 ^(k)
c                            = exp[dm] Delta Lambda(k) Lambda_n

      call quarot( tl(5) , u1 )

c     Copy rotation at t_n+1

      call pmove ( tl    , t1    , 4  )

c     DYNAMIC UPDATES FOR ITERATION WITHIN TIME STEP

      if (fl(9)) then

c       Intermediate Velocity and Acceleration Vectors:

        do i = 1,3
          ws    =   vn(i) + 0.5d0*dt*an(i)
          vs(i) =  -vn(i) + (2.d0 - gamma/beta) * ws
          as(i) =   an(i) - 1.d0/(beta*dt)      * ws
        end do ! i

c       Transform Intermediate Vectors:

        call quavec ( tl(5) , vs , vr )
        call quavec ( tl(5) , as , ar )

c       Update Velocities and Accelerations:

        facv     = gamma / (beta *dt)
        faca     = 1.d0  / (beta *dt*dt)
        do i = 1,3
          v1(i) = facv * u1(i) + vr(i)
          a1(i) = faca * u1(i) + ar(i)
        end do ! i

      endif

      end
