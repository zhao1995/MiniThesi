c$Id:$
      subroutine dyna10(du,urate,nneq,ndf,ndp,ndo,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Move backup storage to 'nrt' based index         20/10/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform static and ODE updates using Euler implicit
c               method for 2nd order problems

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
      include  'part0.h'
      include  'tdata.h'

      integer   i, n, nneq,ndf,isw, ndp(*),ndo(*)
      real*8    cn1,cn2, du(*),urate(nneq,*)

      save

c     Update solution vectors to begin a step

c     Conserving GEMM update
c          u(n,1)     = d(t-n+1);
c          u(n,2)     = d(t-n+1) - d(t-n);
c          u(n,3)     = du(t-n+1)
c                     = [primary dependent variable for solver]
c          urate(n,1) = v(t-n+1);
c          urate(n,2) = a(t-n+1)

c     Starting condition based on d(t-n+1) = d(t-n)

      cn1 = 1.d0/dt
      if(isw.eq.1) then

        do i = 1,ndf
          if(ndp(i).eq.npart .and. ndo(i).ge.1) then
            do n = i,nneq,ndf
              urate(n,nrt-1) =  urate(n,1)
              urate(n,nrt  ) =  urate(n,2)

              urate(n,2)     =  -cn1*urate(n,1)
              urate(n,1)     =   0.0d0
            end do ! n
          endif
        end do ! i

c     Update within step

      elseif(isw.eq.2) then

        cn2 = cn1*cn1
        do i = 1,ndf
          if(ndp(i).eq.npart) then
            do n = i,nneq,ndf
              urate(n,1) = urate(n,1) + cn1*du(n)
              urate(n,2) = urate(n,2) + cn2*du(n)
            end do ! n
          endif
        end do ! i

c     Back up step

      elseif(isw.eq.3) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            do n = i,nneq,ndf
              urate(n,1) =  urate(n,nrt-1)
              urate(n,2) =  urate(n,nrt  )
            end do ! n
          endif
        end do ! i

      endif

      end
