c$Id:$
      subroutine dyna04(du,urate,nneq,ndf,ndp,ndo,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Move backup storage to 'nrt' based index         20/10/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform 2nd order ODE update using explicit Newmark
c               method.

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
c         urate(nneq,nn)    Rate vectors fixed by ALGO
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'comblk.h'
      include  'idptr.h'
      include  'part0.h'
      include  'pointer.h'
      include  'tdata.h'

      integer   i, n, nneq,nneq2,ndf,isw, ndp(*),ndo(*)
      real*8    ub, du(*),urate(nneq,*)

      save

c     Update solution vectors at start of step

      if(isw.eq.1) then

        nneq2 = nneq + nneq
        do i = 1,ndf
          if(ndp(i).eq.npart .and. ndo(i).ge.2) then
            do n = i,nneq,ndf

c             Compute values from forced inputs for fixed dof
c                  ID = np(31)
              if(mr(id31+n-1).gt.0) then

                ub          = dt*urate(n,1) + c1*urate(n,2)
                du(n)       = du(n)         + ub
                du(n+nneq)  =                 ub
                du(n+nneq2) =                 ub

c                 Save velocity and acceleration from t_n

                urate(n,nrt-1) = urate(n,1)
                urate(n,nrt  ) = urate(n,2)

                urate(n,1)     = urate(n,1)    + (dt - c2)*urate(n,2)
                urate(n,2)     = 0.0d0

c             Compute values from current solution state

              else
c                            FTN = np(30)
                du(n)       = hr(np(30)+n-1)
                du(n+nneq)  = du(n) - hr(np(30)+n+nneq-1)
                du(n+nneq2) = du(n+nneq)

c       WARNING: Boundary velocity and accelerations are NOT computed.
c                Do not use consistent mass or problems with damping
c                with specified boundary motions!

              endif
            end do ! n
          endif
        end do ! i

c     Update solution vectors within step

      elseif(isw.eq.2) then

        do i = 1,ndf
          if(ndp(i).eq.npart .and. ndo(i).ge.2) then
            do n = i,nneq,ndf
              urate(n,1) = urate(n,1) + c2*du(n)
              urate(n,2) = urate(n,2) +    du(n)
              du(n)      = 0.0d0
            end do ! n
          endif
        end do ! i

c     Backup solution vectors

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
