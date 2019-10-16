c$Id:$
      subroutine dyna01(du,urate,nneq,ndf,ndp,ndo,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add start options based on 'idyn0' in tdata.h    02/02/2007
c          idyn0 = 1: u_n+1^0 = u_n; 2: v_n+1^0 = v_n; 3: a_n+1^0 = 0
c       2. Change idyn0 = 3 to a_n+1^0 = 0                  04/04/2007
c       3. Move backup storage to 'nrt' based index         20/10/2010
c       4. Add idyn0 = 4: a_n+1^0 = a_n                     14/05/2012
c       5. Remove initializations of 'du' & nnq2            18/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform static and ODE updates using Newmark formulas.

c      Inputs:
c         du(*)             Displacement                    (isw = 1)
c                           Last displacement increment     (isw = 2)
c                           U_n+1 - U_n                     (isw = 3)
c         urate(nneq,*)     Rate vectors - fixed by ALGO
c         nneq              numnp * ndf
c         ndf               Number of DOF/node
c         ndp(*)            Partition dof's
c         ndo(*)            Order dof's: 0 = static, n = n-th order.
c         isw               Control switch
c                            1  STARTING update: begining of time step
c                            2  UPDATE at an iteration within time step
c                            3  BACK solution to begining of time step

c      Outputs:
c         urate(nneq,nn)    Rate vectors fixed by ALGO
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'part0.h'
      include  'ddata.h'
      include  'tdata.h'

      integer   i, n, nneq,ndf,isw, ndp(*),ndo(*)
      real*8    ur1,c6, du(*),urate(nneq,*)

      save

c     Newmark-Beta updates: urate(n,1) = velocity
c                           urate(n,2) = acceleration

c     (1) Update solution vectors to begin a step

      if(isw.eq.1) then

c       Initialize step

        c6   = dt*c1
        do i = 1,ndf
          if(ndp(i).eq.npart .and. ndo(i).ge.1) then
            do n = i,nneq,ndf
              urate(n,nrt-1) = urate(n,1)
              urate(n,nrt  ) = urate(n,2)
            end do ! n
            if(idyn0.eq.1) then              ! u_n+1^0 = u_n
              do n = i,nneq,ndf
                ur1        = - c6*urate(n,1) + c3*urate(n,2)
                urate(n,1) =   c4*urate(n,1) + c5*urate(n,2)
                urate(n,2) =   ur1
              end do ! n
            elseif(idyn0.eq.2) then          ! v_n+1^0 = v_n
              do n = i,nneq,ndf
                du(n+nneq) =  dt*(urate(n,1)
     &                           + (0.5d0 - theta(1))*dt*urate(n,2))
                urate(n,2) =  (theta(2) - 1.d0)/theta(2)*urate(n,2)
                du(n+nneq) =  du(n+nneq) + theta(2)*dt*dt*urate(n,2)
                du(n     ) =  du(n+nneq) + du(n)
              end do ! n
            elseif(idyn0.eq.3) then          ! a_n+1^0 = 0
              ur1 = (0.5d0 - theta(1))*dt
              c6  = (1.0d0 - theta(2))*dt
              do n = i,nneq,ndf
                du(n+nneq) = dt*(urate(n,1) + ur1*urate(n,2))
                du(n     ) = du(n+nneq) + du(n)
                urate(n,1) = urate(n,1) + c6 *urate(n,2)
                urate(n,2) = 0.0d0
              end do ! n
            elseif(idyn0.eq.4) then          ! a_n+1^0 = a_n
              ur1 = 0.5d0*dt
              do n = i,nneq,ndf
                du(n+nneq) = dt*(urate(n,1) + ur1*urate(n,2))
                du(n     ) = du(n+nneq) + du(n)
                urate(n,1) = urate(n,1) + dt*urate(n,2)
              end do ! n
            endif
          endif
        end do ! i

c     (2) Update with in time step

      elseif(isw.eq.2) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            if(ndo(i).ge.2) then
              do n = i,nneq,ndf
                urate(n,1) = urate(n,1) + c2*du(n)
                urate(n,2) = urate(n,2) + c1*du(n)
              end do ! n
            endif
            if(ndo(i).eq.1) then
              do n = i,nneq,ndf
                urate(n,1) = urate(n,1) + c2*du(n)
                urate(n,2) = 0.0d0
              end do ! n
            endif
          endif
        end do ! i

c     (3) Backup solution vectors to reinitiate a step

      elseif(isw.eq.3) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            do n = i,nneq,ndf
              urate(n,1) = urate(n,nrt-1)
              urate(n,2) = urate(n,nrt)
            end do ! n
          endif
        end do ! i

      endif

      end
