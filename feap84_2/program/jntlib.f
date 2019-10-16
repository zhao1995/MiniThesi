c$Id:$
      subroutine jntlib(numjts,jnt,irb,x,u,rcg,rlam,rjt,ebig,
     &                  nrbdof,ndm,ndf,neqs,
     &                  alfl,aufl,bfl,jp,ad,al,au,b)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form arrays for joints connecting rigid bodies.

c               Joint library

c                 Type         Description
c               --------     ---------------
c                   1         Ball & Socket
c                   2         Revolute
c                   3         Slider
c                   4         Plane
c                   5         Translation
c                   6         Angle Control
c                   7         Displacement Control

c      Inputs:
c         numjts     - Number of joints
c         jnt(*)     - Joint description data
c         irb(*,*)   - Equation numbers for rigid body terms
c         x(ndm,*)   - Nodal coordinate array
c         u(ndf,*)   - Current solution at t_n+1
c         rcg(*,*,*) - Array describing center of mass locations & rates
c         rlam(*,*,*)- Array describing rotation array locations & rates
c         rjt(5,*,*) - Array of solution values for lagrange multipliers
c         ebig(3,3,*)- Array of axis orientations for rigid body at t_0
c         nrbdof     - Number of dof/rigid bodies
c         ndm        - Spatial dimension of mesh
c         ndf        - Number dof/node
c         neqs       - Number of symmetric equations in tangent array
c         alfl       - Flag true if joint equations are unsymmetric
c         aufl       - Flag true if joint equations to be assembled
c         bfl        - Flag true if residual to be assembled
c         jp(*)      - Pointer do row/column ends in profile of A

c      Outputs:
c         ad(*)      - Diagonals of tangent array
c         al(*)      - Lower part of tangent array
c         au(*)      - Upper part of tangent array
c         b(*)       - Residual array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'idptr.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   alfl,aufl,bfl
      integer   i,j1,j2, n1, nrev, bnum, n,numjts,nrbdof,ndm,ndf,neqs
      integer   jnt(6,numjts), irb(nrbdof,*), jp(*), ld(3,7)
      real*8    x1(3),x2(3),xh1(3),xh2(3),y1(3),y2(3),yh1(3),yh2(3)
      real*8    r1(3),r2(3),rh1(3),rh2(3),lam(5),pr(3,7), sr(21,21)
      real*8    x(3,3,*),u(ndf,*),rcg(3,11,*)
      real*8    rlam(9,6,*),rjt(5,numjts,*),ebig(3,3,*)
      real*8    ad(*),al(*),au(*),b(*),angle(3), alpr

      save

c     Compute joint forces: not for reactions

      nrev = 0
      do n = 1,numjts

        call pzeroi(ld,21)
        call pzero (pr,21)
        call pzero (sr,441)
        j1 = jnt(2,n)
        j2 = jnt(3,n)

c       Set joint position vectors for first point

c       Rigid body point

        alpr = 1.d0/theta(3)

        if(j1.gt.0) then
          call rxnod(xh1,x1,rh1,r1,yh1,y1,x(1,1,n),rcg,rlam,
     &               alpr,j1,ndm)

c       Fixed point in inertia reference frame

        elseif(j1.eq.0) then
          do i = 1,ndm
            x1(i)  = x(i,1,n)
            y1(i)  = 0.0d0
            yh1(i) = 0.0d0
          end do ! i
        end if

c       Set joint position vectors for second point

        if(jnt(1,n).ne.2.and.jnt(1,n).ne.6) then

          if(jnt(1,n).eq.1 .or. jnt(1,n).eq.4) then
            bnum = 1
          else
            bnum = 2
          endif

c         Rigid body point

          if(j2.gt.0) then
            call rxnod(xh2,x2,rh2,r2,yh2,y2,x(1,bnum,n),rcg,rlam,
     &                 alpr,j2,ndm)

c         Fixed point in inertia reference frame

          elseif(j2.eq.0) then
            do i = 1,ndm
              x2(i)  = x(i,bnum,n)
              y2(i)  = 0.0d0
              yh2(i) = 0.0d0
            end do ! i
          end if

        end if

        alpr = theta(3)

c       Ball & Socket

        if(jnt(1,n).eq.1) then

c         Set joint numbers for first point

          if(j1.gt.0) then
            do i = 1,ndm
              ld(i,1) = irb(i  ,j1)
            end do ! i
            if(ndm.eq.2) then
              ld(3,2) = irb(3,j1)
            elseif(ndm.eq.3) then
              do i = 1,ndm
                ld(i,2) = irb(i+3,j1)
              end do ! i
            endif
          endif

c         Set joint numbers for second point

c         Rigid body point

          if(j2.gt.0) then
            do i = 1,ndm
              ld(i,3) = irb(i  ,j2)
            end do ! i
            if(ndm.eq.2) then
              ld(3,4) = irb(3,j2)
            elseif(ndm.eq.3) then
              do i = 1,ndm
                ld(i,4) = irb(i+3,j2)
              end do ! i
            endif

c         Flexible body point

          elseif(j2.lt.0) then

            n1 = abs(j2) - 1
            call defnod(hr(np(43)+n1*ndm),u(1,-j2),mr(id31+n1*ndf),
     &                   ld(1,3),x2)
            do i = 1,3
              y2(i)  = 0.0d0
              yh2(i) = 0.0d0
            end do ! i

c         Attach to "earth"

          else
            do i = 1,3
              x2(i)  = x(i,1,n)
              y2(i)  = 0.0d0
              yh2(i) = 0.0d0
            end do ! i
          endif

c         Lagrange multiplier quantity

          n1 = jnt(5,n)
          do i = 1,n1
            ld(i,5) = jnt(6,n) + i
            lam(i)  = rjt(i,n,1)
          end do ! i

          if(j1.ne.j2) then
            call jbalsk(lam,x1,x2,y1,y2,yh1,yh2,alpr,pr,sr,alfl)
          else
            write(*,2001) n
          endif

c       Revolute

        elseif(jnt(1,n).eq.2) then

          nrev = nrev + 1

          if(j1.ne.j2) then

c           Rigid body equations

            if(j1.gt.0) then
              do i = 1,3
                ld(i,1) = irb(i+3,j1)
              end do ! i
            end if

            if(j2.gt.0) then
              do i = 1,3
                ld(i,2) = irb(i+3,j2)
              end do ! i
            end if

c           Lagrange multiplier equations

            do i = 1,2
              ld(i,3) = jnt(6,n) + i
              lam(i)  = rjt(i,n,1)
            end do ! i

            call jrevol(lam,rlam,ebig(1,1,nrev),alpr,j1,j2,pr,sr)

          else
            write(*,2001) n
          endif

c       Basic Constraint Type 2

        elseif(jnt(1,n).eq.8) then

          nrev = nrev + 1

          if(j1.ne.j2) then

c           Rigid body equations

            if(j1.gt.0) then
              do i = 1,3
                ld(i,1) = irb(i+3,j1)
              end do ! i
            end if

            if(j2.gt.0) then
              do i = 1,3
                ld(i,2) = irb(i+3,j2)
              end do ! i
            end if
c           Lagrange multiplier equations

            ld(1,3) = jnt(6,n) + 1
            lam(1)  = rjt(1,n,1)

            call jrevo0(lam(1),x(1,1,n),rlam,ebig(1,1,nrev),alpr,
     &                  j1,j2,pr,sr)

          else
            write(*,2001) n
          end if

c       Slider

        elseif(jnt(1,n).eq.3) then

          nrev = nrev + 1

          if(j1.ne.j2) then

c           Set joint numbers for first point

c           Rigid body point

            if(j1.gt.0) then
              do i = 1,3
                ld(i,1) = irb(i  ,j1)
                ld(i,2) = irb(i+3,j1)
              end do ! i
            endif

c           Set joint numbers for second point

c           Rigid body point

            if(j2.gt.0) then
              do i = 1,3
                ld(i,3) = irb(i  ,j2)
                ld(i,4) = irb(i+3,j2)
              end do ! i
            endif

c           Lagrange multiplier equations

            do i = 1,3
              ld(i,5) = jnt(6,n) + i
              lam(i)  = rjt(i,n,1)
            end do ! i
            ld(1,6) = jnt(6,n) + 4
            lam(4)  = rjt(4,n,1)

            call jslide(lam,rlam,ebig(1,1,nrev),alpr,x1,
     &            xh2,x2,rh1,r1,rh2,r2,j1,j2,pr,sr)

          else
            write(*,2001) n
          endif

c       Plane

        elseif(jnt(1,n).eq.4) then

          nrev = nrev + 1

          if(j1.ne.j2) then

c           Rigid body equations

            do i = 1,3
              ld(i,1) = irb(i  ,j1)
              ld(i,2) = irb(i+3,j1)
              ld(i,3) = irb(i  ,j2)
              ld(i,4) = irb(i+3,j2)
            end do ! i

c           Lagrange multiplier equations

            ld(1,5) = jnt(6,n) + 1
            lam(1)  = rjt(1,n,1)

            call jplane(lam,rlam,ebig(1,1,nrev),alpr,x1,
     &            xh2,x2,rh1,r1,rh2,r2,j1,pr,sr)

          else
            write(*,2001) n
          endif

c       Translation

        elseif(jnt(1,n).eq.5) then

          nrev = nrev + 1

          if(j1.ne.j2) then

c           Rigid body equations

            do i = 1,3
              ld(i,1) = irb(i  ,j1)
              ld(i,2) = irb(i+3,j1)
              ld(i,3) = irb(i  ,j2)
              ld(i,4) = irb(i+3,j2)
            end do ! i

c           Lagrange multiplier equations

            do i = 1,5
              ld(i,5) = jnt(6,n) + i
              lam(i)  = rjt(i,n,1)
            end do ! i

            call jtrans(lam,rlam,ebig(1,1,nrev),alpr,x1,
     &            xh2,x2,rh1,r1,rh2,r2,j1,j2,pr,sr)

          else
            write(*,2001) n
          endif

c       Angle Control

        elseif(jnt(1,n).eq.6) then

          nrev = nrev + 1

          if(j1.ne.j2) then

c           Rigid body equations

            if(j1.gt.0) then
              do i = 1,3
                ld(i,1) = irb(i+3,j1)
              end do ! i
            end if

            if(j2.gt.0) then
              do i = 1,3
                ld(i,2) = irb(i+3,j2)
              end do ! i
            end if

            call jangle(j1,j2,rlam,ebig(1,1,nrev),angle)

c           x(3,3,n) = 0  => angle control
c           x(3,3,n) = 1  => torque control
c           x(3,3,n) = 2  => linear elastic spring

            if(x(3,3,n).eq.0) then

c             Lagrange multiplier equations

              ld(1,3) = jnt(6,n) + 1
              lam(1)  = rjt(1,n,1)

              call jangl0(lam(1),angle(2),x(1,3,n),rlam,ebig(1,1,nrev),
     &                    j1,j2,pr,sr)

            elseif(x(3,3,n).eq.1) then

              call jangl1(x(1,3,n),rlam,ebig(1,1,nrev),
     &                    j1,j2,pr,sr)

            elseif(x(3,3,n).eq.2) then

              call jangl2(x(1,3,n),rlam,ebig(1,1,nrev),j1,j2,pr,sr)

            end if

          else
            write(*,2001) n
          end if

c       Displacement Control

        elseif(jnt(1,n).eq.7) then

          nrev = nrev + 1

          if(j1.ne.j2) then

c           Set joint numbers for first point

c           Rigid body point

            if(j1.gt.0) then
              do i = 1,3
                ld(i,1) = irb(i  ,j1)
                ld(i,2) = irb(i+3,j1)
              end do ! i
            end if

c           Set joint numbers for second point

c           Rigid body point

            if(j2.gt.0) then
              do i = 1,3
                ld(i,3) = irb(i  ,j2)
                ld(i,4) = irb(i+3,j2)
              end do ! i
            end if

c           x(3,3,n) = 0  => displacement control
c           x(3,3,n) = 1  => force control
c           x(3,3,n) = 2  => linear elastic spring

            if(x(3,3,n).eq.0) then

c             Lagrange multiplier equations

              ld(1,5) = jnt(6,n) + 1
              lam(1)  = rjt(1,n,1)

              call jdisp0(lam(1),x(1,3,n),rlam,ebig(1,1,nrev),x1,
     &                   xh2,x2,rh1,r1,rh2,r2,j1,pr,sr)

            elseif(x(3,3,n).eq.1) then

              call jdisp1(x(1,3,n),rlam,ebig(1,1,nrev),x1,x2,
     &                    j1,j2,pr,sr)

            elseif(x(3,3,n).eq.2) then

              call jdisp2(x(1,1,n),rlam,ebig(1,1,nrev),
     &                    x1,x2,xh1,xh2,
     &                    r1,r2,rh1,rh2,
     &                    j1,j2,pr,sr)

            end if

          else
            if(ior.lt.0) then
              write(*,2001) n
            endif
            write(ilg,2001) n
            write(iow,2001) n
          end if

        else
          if(ior.lt.0) then
            write(*,2000) n
          endif
          write(ilg,2000) n
          write(iow,2000) n
        end if

        call dasble(sr,pr,ld,jp,21,neqs,aufl,bfl,b,al,au,ad)

      end do ! n

2000  format(' *ERROR* Joint',i3,' not an admissible type')
2001  format(' *ERROR* Joint',i3,' does not connect rigid bodies')

      end
