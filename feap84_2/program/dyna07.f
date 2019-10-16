c$Id:$
      subroutine dyna07(du,urate,nneq,ndf,ndp,ndo,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise location for update of 'urate'            06/07/2007
c       2. Move backup storage to 'nrt' based index         20/10/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform static & 1st order ODE update using generalized
c               midpoint method.

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

      integer   i, n, nneq,ndf,isw

      integer   ndp(*),ndo(*)
      real*8    du(*),urate(nneq,*)

      save

c     FIRST ORDER: Generalized mid-point update

c     Update solution vectors to begin a step

      if(isw.eq.1) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            if(ndo(i).ge.1) then
              do n = i,nneq,ndf
                urate(n,nrt-1) = urate(n,1)
                urate(n,nrt  ) = urate(n,2)
                urate(n,1)     = 0.0d0
                urate(n,2)     = 0.0d0
                urate(n,3)     = du(n)
              end do ! n
            else
              do n = i,nneq,ndf
                urate(n,1) = 0.0d0
                urate(n,2) = 0.0d0
                urate(n,3) = du(n)
              end do ! n
            endif
          endif
        end do ! i

c     Update solution vectors within step

      elseif(isw.eq.2) then

        do i = 1,ndf
          if(ndp(i).eq.npart)then
            if(ndo(i).ge.1) then
              do n = i,nneq,ndf
                urate(n,1) = urate(n,1) + c1*du(n)
                urate(n,2) = urate(n,1)
                urate(n,3) = urate(n,3) + c3*du(n)
              end do ! n
            else
              do n = i,nneq,ndf
                urate(n,3) = urate(n,3) + c3*du(n)
              end do ! n
            endif
          endif
        end do ! i

c     Backup solution vector to reinitiate step

      elseif(isw.eq.3) then

        do i = 1,ndf
          if(ndp(i).eq.npart .and. ndo(i).ge.1) then
            do n = i,nneq,ndf
              urate(n,1) = urate(n,nrt-1)
              urate(n,2) = urate(n,nrt  )
            end do ! n
          endif
        end do ! i

      endif

      end
