c$Id:$
      subroutine dyna02(du,urate,nneq,ndf,ndp,ndo,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Move backup storage to 'nrt' based index         20/10/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Static, 1st and 2nd order ODE integration by backward
c               Euler method.

c      Notes: 1. Values of velocity and acceleration are returned as
c                zero for values of ndo(i) .gt. order specified.

c      Inputs:
c         du(*)             Increment to displacement
c         urate(nneq,*)     Rate vectors - fixed by ALGO
c         nneq              numnp * ndf
c         ndf               Number of DOF/node
c         ndp(*)            Partition dof's
c         ndo(*)            Order dof's
c         isw               Control switch
c                            1  STARTING update: begining of time step
c                            2  UPDATE at an iteration within time step
c                            3  BACK solution to begining of time step

c      Outputs:
c         urate(nneq,nn)    Rate vectors:
c                            1  Velocity at t_n+1     (ndo(i) .ge. 1
c                            2  Acceleration at t_n+1 (ndo(i) .ge. 2
c                            6  Velocity at t_n       (ndo(i) .ge. 1
c                            7  Acceleration at t_n   (ndo(i) .ge. 2
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'part0.h'
      include  'tdata.h'

      integer   i, n, nneq,ndf,isw, ndp(*),ndo(*)
      real*8    du(*),urate(nneq,*)

      save

c     Backward Euler: Initialize at start of step

      if(isw.eq.1) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            do n = i,nneq,ndf
              urate(n,nrt-1) = urate(n,1)
              urate(n,nrt  ) = urate(n,2)
            end do ! n
            if(ndo(i).ge.2) then
              do n = i,nneq,ndf
                urate(n,1) = -   urate(n,1)
                urate(n,2) =     urate(n,1)
              end do ! n
            elseif(ndo(i).eq.1) then
              do n = i,nneq,ndf
                urate(n,1) = 0.0d0
              end do ! n
            endif
          endif
        end do ! i

c     Backward Euler: Updates in iterations

      elseif(isw.eq.2) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            if(ndo(i).ge.2) then
              do n = i,nneq,ndf
                urate(n,1) = urate(n,1) + c1*du(n)
                urate(n,2) = urate(n,1)
              end do ! n
            elseif(ndo(i).eq.1) then
              do n = i,nneq,ndf
                urate(n,1) = urate(n,1) + c1*du(n)
                urate(n,2) = 0.0d0
              end do ! n
            endif
          endif
        end do ! i

c     Backward Euler: Back up to start of step

      elseif(isw.eq.3) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            do n = i,nneq,ndf
              urate(n,1) = urate(n,nrt-1)
              urate(n,2) = urate(n,nrt  )
            end do ! n
          endif
        end do ! i

      endif

      end
