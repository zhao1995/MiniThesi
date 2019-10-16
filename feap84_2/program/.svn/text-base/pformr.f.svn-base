c$Id:$
      subroutine pformr(ld,p,s,x,irb,jp,fn,u,b,a,al,
     &                  aufl,bfl,alfl,dfl,ndm,ndf,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Assemble rigid body arrays into global arrays

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         irb(*,*) -  Rigid body equation numbers
c         jp(*)     - Pointers for row/columns in profile
c         fn(*)     - Nodal force/displacements of mesh
c         u(*)      - Nodal solution values for mesh
c         aufl      - Flag, assemble tangent if true
c         bfl       - Flag, assemble residual if true
c         alfl      - Flag, assemble unsymmetric tangent if true
c         dfl       - Flag, assemble reactions if true
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         nst       - Size of element array
c         isw       - Switch to control quantity to be computed

c      Scratch:
c         ld(*)     - Local/Global equation numbers for element
c         p(*)      - Element residual
c         s(nst,*)  - Element tangent

c      Outputs:
c         b(*)      - Assembled residual for rigid body
c         a(*)      - Assembled diagonal and upper part of tangent
c         al(*)     - Assembled lower part of tangent
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'ddata.h'
      include  'eltran.h'
      include  'eqsym.h'
      include  'evdata.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'modreg.h'
      include  'pointer.h'
      include  'prld1.h'
      include  'prlod.h'
      include  'ptdat6.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'tdata.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   aufl,alfl,bfl,dfl,ldflg, tanfl
      integer   i, ixt, j,k, n, ndm,ndf,nst,isw
      integer   ld(nst),jp(*),irb(nrbdof,*)
      real*8    p(nst),s(nst,nst),x(ndm,*),fn(ndf,numnp,*)
      real*8    u(ndf,numnp,*),b(*),a(*),al(*)
      real*8    rotv(3),rota(3),w1(3),wn(3),pi1(3),pin(3),inert1(3,3)
      real*8    r(3),ra(3),rangm(3), pdotr,thn,alp,rener,propt

      save

c     Rigid body effects: loads, joints, mass effects

c     This version handles the rotational update cases: Pointer = np(96)
c         nreqrb =  3: Cayley transform (THETA_n+a)
c         nreqrb = -3: Cayley transform (THETA_n+1)
c         nreqrb = -4: Exponential map
c         nreqrb = -6: Central Difference explicit

c     Inertia / R arrays
c         RLIST = list of rigid bodies by materials  [       nummat]
c         RIRB  = equation numbers for rigid body    [nrbdof,nummat]
c         RMASS = rigid body mass                    [1,     nummat]
c         RINER = rigid body inertia tensor          [3, 3,  nummat]

c     Order: R_0, r_n, Dr_n, dr_n, rt_n, rtt_n, r_alp, rt_alp, rtt_alp

c     Compute arrays for applied loads, joints, and rigid body mass

      if(isw.eq.3 .or. (.not.dfl .and. isw.eq.6) .or. isw.eq.13) then

        call pzeroi(ld,nst)

c       Loads from fe nodes

        do n = 1,numnp ! {

          ixt = mr(np(100)+n-1)

          if(ixt.gt.0) then

c           Set 'alp' for rotation update method: alp=1/theta: nreqrb 3
c                                                 alp=theta  : nreqrb 5
c                                                 alp=1.0    : nreqrb 7
c                                                 alp=1.0    : nreqrb 8
            if(mr(np(96)+ixt-1).eq.3) then
              alp   = 1.d0/theta(3)
              tanfl = .true.
            elseif(mr(np(96)+ixt-1).eq.-3) then
              alp   = theta(3)
              tanfl = .true.
            elseif(mr(np(96)+ixt-1).eq.-6 .or.
     &             mr(np(96)+ixt-1).eq.-7) then
              alp   = 1.d0
              tanfl = .false.
            else
              alp   = 1.d0
              tanfl = .true.
            endif

c           Look for a nonzero nodal force

            thn   = 1.d0 - theta(3)
            ldflg = .false.
            do i = 1,ndf ! {
              p(i) = theta(3)*fn(i,n,3) + thn*fn(i,n,4)
              if(p(i).ne.0.0d0) ldflg = .true.
            end do ! i     }

c           Compute moment

            if(ldflg) then

              if(tanfl) then
                call rnodld(ra,r,x(1,n),hr(np(95)),ixt,nrk,ndm,ndf,
     &                      hr(np(104)),hr(np(183)),n,mr(np(176)))
              else
                call rnodex(ra,r,x(1,n),hr(np(95)),hr(np(104)),ixt)
              endif

              if(ndm.eq.2) then
                p(3) = p(2)*ra(1) - p(1)*ra(2)
              elseif(ndm.eq.3) then
                p(4) = p(3)*ra(2) - p(2)*ra(3)
                p(5) = p(1)*ra(3) - p(3)*ra(1)
                p(6) = p(2)*ra(1) - p(1)*ra(2)
              endif

c             Compute work done by forces on rigid body

              do i = 1,ndm
                epl(9) = epl(9) + fn(i,n,3)*u(i,n,1)
              end do ! i

c             Set address for assembly

              do i = 1,nrbdof ! {
                ld(i) = irb(i,ixt)
              end do ! i        }

c             Inertial fixed forces from fem nodes

              call pzero(s,nst*nst)
              if(tanfl) then
                if(ndm.eq.2) then
                  s(3,3) = s(3,3) + alp*(p(1)*r(1) + p(2)*r(2))
                elseif(ndm.eq.3) then
                  pdotr = alp*(p(1)*r(1) + p(2)*r(2) + p(3)*r(3))
                  do i = 4,6
                    s(i,i) = s(i,i) + pdotr
                    if(alfl) then
                      do j = 4,6
                        s(i,j) = s(i,j) - alp*r(i-3)*p(j-3)
                      end do ! j
                    else
                      do j = 4,6
                        s(i,j) = s(i,j) - alp*(r(i-3)*p(j-3)
     &                                  +      r(j-3)*p(i-3))*0.5d0
                      end do ! j
                    end if
                  end do ! i
                endif
              endif

c             Assemble rigid body load equations

              call dasble(s,p,ld,jp,nst,neqs,aufl,bfl,
     &                    b,al,a(neq+1),a)
            endif

          endif

        end do ! n       }

c       Assemble direct loads on rigid body (no tangent)

        propt = (1.d0 - theta(3))*propo + theta(3)*prop
        do n = 1,nrlds ! {

          fp(1)= np(106)  + n*7 - 7
          ixt = nint(hr(fp(1)))
          do i = 1,nrbdof ! {
            ld(i) = irb(i,ixt)
            p(i)  = hr(fp(1)+i)*propt
          end do ! i        }

          call dasble(s,p,ld,jp,nrbdof,neqs,.false.,bfl,b,al,a,a)

        end do ! n       }

c       Impose joint restraints

        if(numjts .gt. 0 .and. ctan(1).gt.0.0d0) then
          call jntlib(numjts,mr(np(101)),mr(np(99)),hr(np(102)),u,
     &                hr(np(95)),hr(np(104)),hr(np(103)),hr(np(97)),
     &                nrbdof,ndm,ndf,neqs,
     &                alfl,aufl,bfl,jp,a,al,a(neq+1),b)
        end if

c       Compute rigid body inertia and tangents and assemble

        if( fl(9) .and. ctan(3).gt.0.0d0 ) then

          call pzeroi(ld,nst)
          do n = 0,nrbody-1 ! {

            if(mr(np(96)+n).ne.-6) then

              call inerrb(hr(np(98)+9*n),hr(np(104)+54*n),theta,inert1,
     &                    rotv,rota,pi1,pin,w1,wn)

              if(mr(np(96)+n).eq.3) then

                call tanrb8(hr(np(107)+n),pi1,ctan,dt,inert1,rotv,nst,s)
                alp  = theta(3)

              elseif(mr(np(96)+n).eq.-3) then

                call tanrb8(hr(np(107)+n),pi1,ctan,dt,inert1,rotv,nst,s)
                alp  = 1.d0

              elseif(mr(np(96)+n).eq.-4) then

                call tanrb7(hr(np(107)+n),pi1,theta,dt,inert1,rotv,rota,
     &                      nst, s)
                alp  = 1.d0

              endif

              call resdrb(hr(np(107)+n),pi1,pin,alp,hr(np(95)+33*n),
     &                    dt, p)

            elseif(mr(np(96)+n).eq.-6) then

              call inerbx(hr(np(107)+n),hr(np(95)+33*n),hr(np(104)+54*n)
     &                   ,hr(np(98)+9*n),ctan, nst,s,p,pi1,w1)

            endif

c           Add modal force contributions to rigid body

            if(nmbody.gt.0) then
              do k = 1,nmbody
                if(modbod(k).eq.n+1 .and. mf.gt.0) then
                  if(ndm.eq.2) then
                    p(3)   = p(3)   + br(1)
                    s(3,3) = s(3,3) + ar(1,1)
                  elseif(ndm.eq.3) then
                    do i = 1,3
                      p(i+3) = p(i+3) + br(i)
                      do j = 1,3
                        s(j+3,i+3) = s(j+3,i+3) + ar(j,i)
                      end do ! j
                    end do ! i
                  endif
                endif
              end do ! k
            endif

c           Assemble rigid body equations

            do i = 1,nrbdof ! {
              ld(i) = irb(i,n+1)
            end do ! i        }

            call dasble(s,p,ld,jp,nst,neqs,aufl,bfl,b,al,a(neq+1),a)

c           Compute energy for rigid body

            call enerrb(hr(np(107)+n),pi1,w1,hr(np(95)+33*n),
     &                  ra,rangm,rener)

c           Accumulate the energy and momenta

            do i = 1,3
              epl(i  ) = epl(i  ) + ra(i)
              epl(i+3) = epl(i+3) + rangm(i)
            end do ! i

            epl(7) = epl(7) + rener

          end do ! n    }

        end if

      end if

      end
