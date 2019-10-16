c$Id:$
      subroutine upmodal(u,phi,wf,umode,ixt,imf,rlam,
     &                   ndm,ndf,numnp,nmod,ud)

c     * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Project Modal solutions to nodal solutions
c              assuming x_n+a = r_n+a + lam_n+a(X-R)
c                             + alph*lam_n+1*wf_n+1
c                             + (1-alph)*lam_n*wf_n

c     Inputs:
c       x(ndm,*)       - Nodal coordinates for mesh
c       phi(neqmf,mf)  - Mode shapes
c       wf(mf,3)       - Modal participations
c       umode(neqmf,3) - Modal displacements
c       ixt(*)         - Indicator for nodes attached to rigid body
c       imf(ndf,numnp) - Equation locators for modes
c       rlam           - Rigid body rotation solutions
c       ndm            - Spatial dimension of mesh
c       ndf            - Number dof/node
c       numnp          - Number nodes in mesh

c     Outputs:
c       u(ndf,numnp,*) - Nodal displacements for rigid body nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'evdata.h'
      include  'iofile.h'
      include  'modreg.h'
      include  'pointer.h'
      include  'tdata.h'
      include  'comblk.h'

      integer   ndm,ndf,numnp,i,j,k,n,nmod, ixt(*),imf(ndf,*)
      real*8    u(ndf,numnp,*),umode(neqmf,8),ud(ndf,numnp,*)
      real*8    phi(neqmf,*), wf(mf,*), rlam(9,6,*)
      real*8    dx1(3),dx2(3),ul(6,8),dv1,dva
      real*8    lam1(3,3), lamn(3,3), lama(3,3)

      save

c     Compute modal displacements

      do i = 1,8
        do n = 1,neqmf
          umode(n,i) = 0.d0
          do k = 1,mf
            umode(n,i) = umode(n,i) + phi(n,k)*wf(k,i)
          end do ! k
        end do ! n
      end do ! i

c     Convert rigid body displacements to nodal displacements

      do n = 1,numnp

        if(ixt(n).eq.modbod(nmod)) then

c         Extract displacements of node

          do i = 1,ndm
            if(imf(i,n).gt.0) then
              do j = 1,8
                ul(i,j) = umode(imf(i,n),j)
              end do ! j
            else
              do j = 1,8
                ul(i,j) = 0.d0
              end do ! j
            endif
          end do ! i

c         Compute total and incremental rotational displacement at t_n+1

          k = ixt(n)

          call quavec(rlam(1,3,k),ul(1,1),dx1)
          call quavec(rlam(1,3,k),ul(1,2),dx2)

c         Compute total displacement at t_n+1, t_n+a
c         and incremental displacement

          do i = 1,ndm
            u(i,n,1)    = u(i,n,1)    + dx1(i)
            u(i,n,2)    = u(i,n,2)    + dx2(i)
            ud(i,n,nrk) = ud(i,n,nrk) + theta(3)*dx1(i)
          end do ! i

c         Compute rotation matrix at t_n+1, t_n+a

          call quamat(rlam(1,3,k),lam1)
          call quamat(rlam(1,1,k),lamn)
          do i = 1,ndm
            do j = 1,ndm
              lama(i,j) = theta(3)*lam1(i,j) + (1.d0-theta(3))*lamn(i,j)
            end do ! j
          end do ! i

c         Compute velocities at t_n+1, t_n+a

          do i = 1,ndm
            dv1 = 0.0d0
            dva = 0.0d0
            do j = 1,ndm
              dv1 = dv1 + lam1(i,j)*ul(j,4)
              dva = dva + lama(i,j)*ul(j,7)
            end do ! j
            ud(i,n,1) = ud(i,n,1) + dv1
            ud(i,n,4) = ud(i,n,4) + dva
          end do ! i

        end if

      end do ! n

      end
