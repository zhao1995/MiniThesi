c$Id:$
      subroutine rupdat(du,rcg,drot,rlam,irb,eqrb,ndm,fdyn,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Update rigid body solution states for static and
c               dynamic problems

c      Inputs:
c         du(*)     - Increment to solution
c         irb(*)    - Rigid body equation numbers
c         eqrb(*)   - Type of rigid body solution
c         ndm       - Spatial dimension of mesh
c         fdyn      - Flag, Dynamic problem if true
c         isw       - Switch: = 1 Initialize for new time step
c                             = 2 Increment solutions within a step
c                             = 3 Back-up step
c      Outputs:
c         rcg(3,*)  - Rigid body solutions for translations
c         drot(*)   - Incremental rotations
c         rlam(9,*) - Rotational solutions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'iofile.h'
      include  'rigid1.h'
      include  'tdata.h'
      include  'tdatb.h'

      logical   fdyn
      integer   isw,i,j,n,ndm,ndf,nrb,nrbq
      integer   irb(nrbdof,*),eqrb(*),ndfp(10),ndfo(10)
      real*8    du(*), rcg(3,11,*),rlam(9,6,*), drot(6,*)

      save

      data      ndf / 3 /

      call pzero(drot,nrbody*6)

      nrbq = ndf*nrbody

      do nrb = 1,ndf ! {
        ndfp(nrb) = nrbprt
        ndfo(nrb) = 10
      end do ! nrb }

c     Update solution vectors to begin a step

      if(isw.eq.1 .and. fdyn) then

        do nrb = 1,nrbody ! {
          do j = 1,3
            rcg(j,3,nrb) = 0.0d0
          end do ! j
        end do        ! nrb }

c         Newmark update

        do nrb = 1, nrbody ! {

          if(noi.eq.1) then

            call dyna01(rcg(1,2,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,1)

c         Backward Euler update

          elseif(noi.eq.2) then

            call dyna02(rcg(1,2,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,1)

c         Conserving HHT update

          elseif(noi.eq.3) then

            call dyna03(rcg(1,2,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,1)

c         Newmark explicit update

          elseif(noi.eq.4) then


c         Three-parameter algorithm in Conservation form

          elseif(noi.eq.5) then

c           rcg(i,1,n) = r(0)
c           rcg(i,2,n) = r(t-n+1)
c           rcg(i,3,n) = r(t-n+a) - r(t-n);
c           rcg(i,4,n) = dd(t_n+1)
c                      = [primary dependent variable for solver]
c           rcg(i,5,n) = v(t-n+1)
c           rcg(i,6,n) = a(t-n+1)
c           rcg(i,7,n) = u(t-n+a)
c           rcg(i,8,n) = v(t-n+a);
c           rcg(i,9,n) = a(t-n+a)
c           rcg(i,10,n)= v(t-n)
c           rcg(i,11,n)= a(t-n)

c           Starting condition based on d(t-n+1) = d(t-n)

c           Save velocity at t_n

            do j = 1,3
              rcg(j,10,nrb) = rcg(j,5,nrb)
            end do ! j

            call dyna05(rcg(1,2,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,1)

c         Static generalized mid-point integrationi

          elseif(noi.eq.6) then

            call dyna06(rcg(1,2,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,1)

c         FIRST ORDER: generalized mid-point integration

          elseif(noi.eq.7) then

            call dyna07(rcg(1,2,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,1)

c         Central Difference: Explicit

          elseif(noi.eq.8) then

            call dyna08(rcg(1,2,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,1)

          endif

        end do        ! nrb }

c     Update for iterations within step

      elseif(isw.eq.2) then

        if(ndm.eq.2) then
          i = 3
        elseif(ndm.eq.3) then
          i = 0
        endif

        do nrb = 1,nrbody
          do n = 1,ndm
            j = irb(n,nrb)
            if (j.gt.0) then

c             For active degrees-of-freedom compute values from solution
c             where 'du(j)' is an increment of 'u' for active dof 'j'.

              rcg(n,2,nrb) = rcg(n,2,nrb) + cc1*du(j)
              rcg(n,3,nrb) = rcg(n,3,nrb) + cc2*du(j)
              rcg(n,4,nrb) =                    du(j)
              if(noi.eq.0) then
                rcg(n,7,nrb) = rcg(n,7,nrb) + cc1*du(j)
              endif
c           else

c             For fixed dofs compute values from forced inputs

c             ub         = 0.0d0

c             Incremental boundary solutions anot needed for updates

c             if(cc1.ne.0.0d0) then
c               rcg(n,4,nrb) = (ub - rcg(n,2,nrb))/cc1
c             else
c               rcg(n,4,nrb) = 0.0d0
c               if(ub-rcg(n,2,nrb).ne.0.0d0) then
c                 write(iow,*)' WARNING - infinite acceleration'
c               endif
c             endif
c             rcg(n,3,nrb)  = rcg(n,3,nrb) + rcg(n,4,nrb)
c             rcg(n,2,nrb)  = rcg(n,1,nrb) + ub

            endif
          end do ! n

c         Store rigid body increments for body

          do n = 1, nrbdof
            j = irb(n,nrb)
            if (j.gt.0) then
              drot(n+i,nrb) = du(j)
            end if
          end do ! n

        end do ! nrb

c       Update transient solutions

        if(fdyn) then

          do nrb = 1,nrbody

c           Update for Newmark

            if (noi.eq.1) then

              call dyna01(rcg(1,4,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,2)

c           Update for Backward Euler

            elseif (noi.eq.2) then

              call dyna02(rcg(1,4,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,2)

c           Update for HHT

            elseif (noi.eq.3) then

              call dyna03(rcg(1,4,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,2)

c           Update for explicit algorithm

            elseif (noi.eq.4) then

c           Update conserving algorithm solution

            elseif(noi.eq.5) then

              call dyna05(rcg(1,4,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,2)

c           Update for STATIC generalized mid-point rule

            elseif(noi.eq.6) then

              call dyna06(rcg(1,4,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,2)

c           Update for FIRST ORDER generalized mid-point rule

            elseif(noi.eq.7) then

              call dyna07(rcg(1,4,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,2)

c           Central Difference: Explicit

            elseif(noi.eq.8) then

              call dyna08(rcg(1,4,nrb),rcg(1,5,nrb),ndf,ndf,ndfp,ndfo,2)

            endif

          end do ! nrb

        endif

c     Backup solution vectors to reinitiate a step

      elseif(isw.eq.3) then

      endif

c     Update rotational degrees of freedom on rigid bodies

      call updrot( drot, 6, rlam, eqrb , nrbody, isw )

      end
