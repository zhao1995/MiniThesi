c$Id:$
      subroutine dyna06(du,urate,nneq,ndf,ndp,ndo,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Move backup storage to 'nrt' based index         20/10/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform static update using generalized midpoint rule.

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

      include  'part0.h'
      include  'tdata.h'

      integer   i, n, nneq,ndf,isw, ndp(*),ndo(*)
      real*8    du(*),urate(nneq,*)

      save

c     STATICS:  Generalized Mid-point configuration

c          u(n,1)     = d(t-n+1)
c          u(n,2)     = d(t-n+a) - d(t-n);
c          u(n,3)     = dd(t-n+1)
c                     = [primary dependent variable for solver]
c          urate(n,3) = u(t-n+a)

c     Starting condition based on d(t-n+1) = d(t-n)

c     Update solution vectors to begin a step

      if(isw.eq.1) then

        do i = 1,ndf
          if(ndp(i).eq.npart .and. ndo(i).ge.0) then
            do n = i,nneq,ndf
              urate(n,3) =  du(n)
            end do ! n
          endif
        end do ! i

c     Conserving algorithm solution

      elseif(isw.eq.2) then

        do i = 1,ndf
          if(ndp(i).eq.npart .and. ndo(i).ge.0) then
            do n = i,nneq,ndf
              urate(n,3) = urate(n,3) + c3*du(n)
            end do ! n
          endif
        end do ! i

c     Backup solution vectors to reinitiate a step

      elseif(isw.eq.3) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            do n = i,nneq,ndf
              urate(n,3) = urate(n,3) - du(n)
            end do ! n
          endif
        end do ! i

      endif

      end
