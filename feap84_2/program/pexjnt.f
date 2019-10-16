c$Id:$
      subroutine pexjnt(jnt,xjnt,ebig, irbd, exin,dr,irb)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Transform joints for explicit solution

c      Inputs:
c         jnt(6,*)    - Joint markers
c         xjnt(3,3,*) - Joint coordinates
c         ebig(3,3,*) - Revolute orientations
c         exin(6,6,*) - Rigid body inertia properties with flexible part

c      Outputs:
c         irbd(*)     - Rigid numbers
c         dr(*)       - Residual
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   n,n1,n2, j,j1,j2,jtype
      integer   jnt(6,numjts), irbd(nrbody), irb(nrbdof,nrbody)
      real*8    xjnt(3,3,numjts),ebig(3,3,numjts)
      real*8    exin(6,6,nrbody),dr(*), resd(6,2),accel(6,2)

      save

      do n = 1,numjts

        jtype    = jnt(1,n)
        n1       = jnt(2,n)
        n2       = jnt(3,n)
        irbd(n1) = 1
        irbd(n2) = 1
        fp(1)     = np( 95) + 33*n - 33
        fp(2)     = np(104) + 54*n - 54
        do j = 1,nrbdof
          j1 = irb(j,n1)
          if(j1.gt.0) then
            resd(j,1) = dr(j1)
          else
            resd(j,1) = 0.0d0
          endif
          j2 = irb(j,n2)
          if(j2.gt.0) then
            resd(j,2) = dr(j2)
          else
            resd(j,2) = 0.0d0
          endif
        end do ! j
        if(jtype.eq.1) then  ! Spherical Joint
          call rsphere(n1,n2, xjnt(1,1,n),hr(fp(1)),hr(fp(2)),exin,resd,
     &                 accel)
        elseif(jtype.eq.2) then  ! Revolute Joint
          call rrevolu(n1,n2, xjnt(1,1,n),hr(fp(1)),hr(fp(2)),
     &                 ebig(1,1,n),exin,resd,accel)
        else
          write(iow,3000)
          call plstop()
        endif
        do j = 1,nrbdof
          j1 = irb(j,n1)
          if(j1.gt.0) then
            dr(j1) = accel(j,1)
          endif
          j2 = irb(j,n2)
          if(j2.gt.0) then
            dr(j2) = accel(j,2)
          endif
        end do ! j
      end do ! n

c     Formats

3000  format(/'  *ERROR* JOINT TYPE NOT IMPLEMENTED IN EXPLICIT')

      end
