c$Id:$
      subroutine rasble(x,u,mass,dr,ixt,irb,rcg,exin,
     &                  jnt,xjnt,ebig)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Transform nodal mass to rigid body inertias

c      Inputs:
c         x(ndm,*)    - Nodal coordinates
c         u(ndf,*)    - Solution vector
c         mass(ndf,*) - Nodal Mass array
c         dr(*)       - Residual
c         id(ndf,*)   - Global equation numbers
c         ixt(*)      - Nodal rigid body numbers (0 = flexible node)
c         irb(*)      - Rigid body equation numbers
c         rcg(3,*,*)  - Rigid body translational solution/rates
c         jnt(6,*)    - Joint markers
c         xjnt(3,3,*) - Joint coordinates
c         ebig(3,3,*) - Revolute orientations

c      Outputs:
c         exin(6,6,*) - Rigid body inertia properties with flexible part
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'gltran.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'sdata.h'
      include  'comblk.h'

      integer   i,j,k,n, ixt(*),irb(nrbdof,*),jnt(6,numjts),irbd(10)
      real*8    mxmass, x(ndm,*), u(ndf,*), mass(ndf,*), dr(*)
      real*8    xjnt(3,3,numjts), ebig(3,3,numjts)
      real*8    rcg(3,11,*), r(6), h1(3,3), exin(6,6,*), du(6)

      save

c     Modify inertia matrix for flexible contribution to rigid body

      do n = 0,nrbody-1
        call inerbx(hr(np(107)+n),hr(np(95)+33*n),hr(np(104)+54*n),
     &              hr(np(98)+9*n),gtan, 6, exin(1,1,n+1),r, du,h1)
        irbd(n+1) = 0
      end do ! n

      do i = 1,numnp

c       Check if node rigid [ ixt(i) > 0 ]

        n = ixt(i)

        if( n.gt.0 ) then

c         Check for mass on interfaces

          mxmass = 0.d0
          do j = 1,ndf
            mxmass = max(mxmass,abs(mass(j,i)))
          end do ! j

          if(mxmass.gt.0.0d0) then

c           Set up relative position vector for rigid node

            r(3)  = 0.0d0
            do j = 1,ndm
              r(j)  = x(j,i) + u(j,i) - rcg(j,2,n)
            end do ! j

c           Translational dof coupling to translational rigid nodes

            do j = 1,ndf
              exin(j,j,n) = exin(j,j,n) + mass(j,i)
            end do ! j

c           Rotational dof coupling to translational rigid nodes

            if(ndm.eq.2) then
              h1(1,1)   = -r(2)*mass(1,i)
              h1(2,1)   =  r(1)*mass(2,i)
              exin(1,3,n) = exin(1,3,n) + h1(1,1)
              exin(2,3,n) = exin(2,3,n) + h1(2,1)
              exin(3,1,n) = exin(3,1,n) + h1(1,1)
              exin(3,2,n) = exin(3,2,n) + h1(2,1)
              exin(3,3,n) = exin(3,3,n) - h1(1,1)*r(2) + h1(2,1)*r(1)
            elseif(ndm.eq.3) then
               h1(1,1)   =  0.0d0
               h1(1,2)   = -r(3)*mass(2,i)
               h1(1,3)   =  r(2)*mass(3,i)
               h1(2,1)   =  r(3)*mass(1,i)
               h1(2,2)   =  0.0d0
               h1(2,3)   = -r(1)*mass(3,i)
               h1(3,1)   = -r(2)*mass(1,i)
               h1(3,2)   =  r(1)*mass(2,i)
               h1(3,3)   =  0.0d0
               do k = 1,3
                 do j = 1,3
                   exin(j+3,k,n) = exin(j+3,k,n) + h1(j,k)
                   exin(k,j+3,n) = exin(k,j+3,n) + h1(j,k)
                 end do ! j
                 exin(k+3,4,n) = exin(k+3,4,n) - h1(k,2)*r(3)
     &                                         + h1(k,3)*r(2)
                 exin(k+3,5,n) = exin(k+3,5,n) - h1(k,3)*r(1)
     &                                         + h1(k,1)*r(3)
                 exin(k+3,6,n) = exin(k+3,6,n) - h1(k,1)*r(2)
     &                                         + h1(k,2)*r(1)
               end do ! k
             end if
          endif
        endif

      end do ! i

c     Solve rigid body parts associated with joints

      call pexjnt(jnt,xjnt,ebig, irbd, exin,dr,irb)

c     Solve rigid body parts not associated with joints

      do n = 1,nrbody
        if(irbd(n).eq.0) then
          call invert (exin(1,1,n),nrbdof,6)   ! invert inertia

          do i = 1,nrbdof
            r(i)  = dr(irb(i,n))               ! load residual
            du(i) = 0.0d0
          end do ! i
          do i = 1,nrbdof
            do j = 1,nrbdof
              du(i) = du(i) + exin(i,j,n)*r(j) ! compute increment
            end do ! j
          end do ! i
          do i = 1,nrbdof
            dr(irb(i,n)) = du(i)               ! move increment
          end do ! i
        endif
      end do ! n

      end
