c$Id:$
      subroutine uprigb(u,urate,x,ixt,rcg,rlam,drot,eqrb,
     &                  ndm,ndf,numnp,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Project rigid body solutions to nodal solutions

c      Inputs:
c         x(ndm,*)       - Nodal coordinates for mesh
c         ixt(*)         - Indicator for nodes attached to rigid body
c         rcg(3,11,*)    - Rigid body translations solutions
c         rlam(3,3,6,*)  - Rigid body rotation solutions
c         drot(*)        - Incremental rotation vector
c         eqrb(*)        - Equations for Rigid Body
c         ndm            - Spatial dimension of mesh
c         ndf            - Number dof/node
c         numnp          - Number nodes in mesh

c      Outputs:
c         u(ndf,numnp,*) - Nodal displacements for rigid body nodes
c         urate(ndf,numnp,*) - Nodal rates for rigid body nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'tdata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   ndm,ndf,numnp,isw,i,j,k,n, ixt(*),eqrb(*)
      real*8    thetn,thet1,du1, rcg(3,11,*), rlam(3,3,6,*), drot(*)
      real*8    u(ndf,numnp,*),urate(ndf,numnp,*),x(ndm,*)
      real*8    dxr(3),dxn(3),dx1(3),dxa(3),omega(3),alpha(3)

      save

c     Convert rigid body displacements to nodal displacements

      thetn  = 1.d0 - theta(3)
      thet1  = theta(3)
      dxr(3) = 0.0d0

      do n = 1,numnp

        if(ixt(n) .gt. 0 ) then

          k = ixt(n)

c         Compute rotational displacement at t_n, t_n+1, and t_n+a

          do i = 1,ndm
            dxr(i) = x(i,n) - rcg(i,1,k)
          end do ! i

c         Multiply rotation matrix by dxr

          if(eqrb(k).eq.-6) then

            do i = 1,3
              dx1(i) = 0.0d0
              do j = 1,3
                dx1(i)   = dx1(i) + rlam(i,j,3,k)*dxr(j)
              end do ! j
              omega(i) = rlam(i,2,5,k)
              alpha(i) = rlam(i,3,5,k)
            end do ! i
            call sasbtdc(omega,omega,dx1, dxn)

            if(abs(isw).eq.1) then
              do i = 1,3
                u(i,n,1) = (rcg(i,2,k) + dx1(i)) - x(i,n)
              end do ! i

              urate(1,n,1) = rcg(1,5,k) + omega(2)*dx1(3)
     &                                  - omega(3)*dx1(2)
              urate(2,n,1) = rcg(2,5,k) + omega(3)*dx1(1)
     &                                  - omega(1)*dx1(3)
              urate(3,n,1) = rcg(3,5,k) + omega(1)*dx1(2)
     &                                  - omega(2)*dx1(1)

              urate(1,n,2) = rcg(1,6,k) + alpha(2)*dx1(3)
     &                                  - alpha(3)*dx1(2)
              urate(2,n,2) = rcg(2,6,k) + alpha(3)*dx1(1)
     &                                  - alpha(1)*dx1(3)
              urate(3,n,2) = rcg(3,6,k) + alpha(1)*dx1(2)
     &                                  - alpha(2)*dx1(1)
              do i = 1,3
                urate(i,n,2) = urate(i,n,2) - dxn(i)
              end do ! i

c           Increment to acceleration at interface node

            elseif(abs(isw).eq.2) then

              u(1,n,3) = rcg(1,6,k) + alpha(2)*dx1(3) - alpha(3)*dx1(2)
              u(2,n,3) = rcg(2,6,k) + alpha(3)*dx1(1) - alpha(1)*dx1(3)
              u(3,n,3) = rcg(3,6,k) + alpha(1)*dx1(2) - alpha(2)*dx1(1)

              do i = 1,3
                u(i,n,3) = u(i,n,3) - dxn(i) - urate(i,n,2)
              end do ! i
            endif

          elseif(abs(isw).eq.2) then

            if(eqrb(k).eq.-7) then
              dx1(1) = dxr(3)*drot(5) - dxr(1)*drot(6)
              dx1(2) = dxr(2)*drot(6) - dxr(3)*drot(4)
              dx1(3) = dxr(1)*drot(4) - dxr(2)*drot(5)
              do i = 1,ndm
                u(i,n,1) = u(i,n,1) + drot(i) + dx1(i)
                u(i,n,2) = u(i,n,2) + drot(i) + dx1(i)
                u(i,n,3) = u(i,n,1) + drot(i) + dx1(i)
              end do ! i
              do i = ndm+1,min(ndf,ndm+3)
                u(i,n,1) = u(i,n,1) + drot(i)
                u(i,n,2) = u(i,n,2) + drot(i)
                u(i,n,3) = u(i,n,1) + drot(i)
              end do ! i

            else
              if(eqrb(k).eq.-6) then
                do i = 1,3
                  dxn(i) = rlam(1,1,1,k)*dxr(1)
     &                   + rlam(2,1,1,k)*dxr(2)
     &                   + rlam(3,1,1,k)*dxr(3)
                  dx1(i) = rlam(1,1,3,k)*dxr(1)
     &                   + rlam(2,1,3,k)*dxr(2)
     &                   + rlam(3,1,3,k)*dxr(3)
                end do ! i
              else
                call quavec(rlam(1,1,1,k),dxr,dxn)
                call quavec(rlam(1,1,3,k),dxr,dx1)
              endif

              do i = 1,ndm
                du1      = u(i,n,1)
                u(i,n,1) = rcg(i,2,k) - x(i,n) + dx1(i)
                u(i,n,3) = u(i,n,1) - du1
              end do ! i
            endif

            if(eqrb(k).eq.3) then

              do i = 1,3
                dxa(i) = thetn * dxn(i) + thet1 * dx1(i)
              end do ! i

            elseif(eqrb(k).eq.-3) then

              do i = 1,3
                dxa(i) = thetn * dxn(i) + thet1 * dx1(i)
              end do ! i

            else

              call quavec(rlam(1,1,2,k),dxr,dxa)

            endif

            if(eqrb(k).ne.-7) then
              do i = 1,ndm
                u(i,n,2) = rcg(i,3,k) + dxa(i) - dxn(i)
              end do ! i
            endif

          endif ! update type

        endif ! rigid body node

      end do ! n

      end
