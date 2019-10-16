c$Id:$
      subroutine frame3d(d,ul,xl,ix,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add implicit/explicit integration option         04/04/2007
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Three dimensional frame element

c     Control data:
c         ndm - 3 (x,y,z)
c         ndf - 6 (u,v,w, theta_x,theta_y,theta_z)
c         nen - 3 or more (see below)

c      Beam end nodes 1 and 2
c      Plane defined by nodes 1, 2, 3 contains z-axis
c                            (perpendicular to x-axis)

c      Vector products: e_1 =  (x_2 - x_1)/|x_2 - x_1|
c                       v_2 = - e_1  x ( x_3 - x_1)
c                       e_2 =   v_2/|v_2|
c                       e_3 =   e_1 x e_2

c                       z (e_3)  x (e_1)
c                     3 o- - - /
c                       |     o 2
c                       |    /
c                       |   / <--- Frame axis
c                       |  /
c                       | /
c                       |/
c     (e_2) y ----------o 1

c     Displacement:     u_x = u_0 + z * theta_y - y * theta_z

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'pmod2d.h'
      include  'comblk.h'

      integer   ndf,ndm,nst,isw, tdof, i, ix(*)
      real*8    xs,ctan1
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),s(nst,*),r(*), pt(3,4)

      save
c     Output descriptor

      if(isw.eq.0) then

        if(ior.lt.0) then
          write(*,*) '   Frame3d: 3-d Frame Element'
        endif

c     Critical time step computation

      elseif(isw.eq.21) then

        call tcritnd(d,xl,ul,ndm,ndf,nel)

c     Compute body force

      elseif(isw.eq.15 .or. isw.eq.23 ) then

        call fbody3d(d,xl, r, ndm,ndf, isw)

c     Check on external nodes

      elseif(isw.eq.26) then

c     Plot local basis triad

      elseif(isw.eq.30) then

        if(nel.eq.2) then

          call beamtri(xl,ndm,d,xs,pt)

          do i = 1,3
            pt(i,2) = pt(i,2)*xs + pt(i,1)
            pt(i,3) = pt(i,3)*xs + pt(i,1)
            pt(i,4) = pt(i,4)*xs + pt(i,1)
          end do !i

c         Plot red for normal 1

          call pppcol(2,1)
          call plotl(pt(1,1),pt(2,1),pt(3,1),3)
          call plotl(pt(1,2),pt(2,2),pt(3,2),2)

c         Plot green for normal 2

          call pppcol(3,1)
          call plotl(pt(1,1),pt(2,1),pt(3,1),3)
          call plotl(pt(1,3),pt(2,3),pt(3,3),2)

c         Plot white for tangent

          call pppcol(1,1)
          call plotl(pt(1,1),pt(2,1),pt(3,1),3)
          call plotl(pt(1,4),pt(2,4),pt(3,4),2)

        endif

c     Call small or finite deformation routine

      else

        if(isw.eq.1) then
          write(iow,2000)
          if(ior.lt.0) write(*,2000)
          call inmate(d,tdof,0,3)

c         Plot line for element

          pstyp = 1

          if(ndf.lt.6.or.ndm.ne.3) then
            write(iow,2001) ndf,ndm
            call plstop()
          end if

c         Deactivate dof in element for dof > 6

          do i = 7,ndf
            ix(i) = 0
          end do ! i

c         Set rotation dof

          do i = 1,2
            ea(i,-iel) = i
            er(i,-iel) = i+3
          end do ! i
          do i = 1,3
            ea3(i,-iel) = i
            er3(i,-iel) = i+3
          end do ! i

c         Check for reference node definition

          if(int(d(96)).lt.1 .or. int(d(96)).gt.4) then
            write(iow,2002) int(d(96))
            call plstop()
          endif

c         Check for inelastic model and no section data

          if(nint(d(40)).gt.0 .and. nint(d(100)).eq.0) then
            write(iow,2003)
            call plstop()
          endif

        endif

c       Explicit/Implicit element solutions

        if(isw.eq.3) then
          ctan1 = ctan(1)
          if(d(187).gt.0.0d0 .and.
     &       min(ctan(2),ctan(3)).gt.0.0d0) then
            ctan(1) = 0.0d0
          endif
        endif

c       Small deformation

        if(d(18).gt.0.0d0) then

c         Shear deformable

          if(nint(d(79)).eq.0) then

c           Shear: 2-node, linear displacement, constant strains

            call frans3l(d,ul,xl,s,r,ndf,ndm,nst,isw)

c         No shear: 2-node, cubic displacement, linear strain

          else
            if(nint(d(100)).eq.0) then
              call frams3d(d,ul,xl,s,r,ndf,ndm,nst,isw)
            else
              call franf3d(d,ul,xl,s,r,ndf,ndm,nst,isw)
            endif
          endif

c       Finite deformation

        else

c         With shear: 2-node, linear, finite strain

          if(nint(d(79)).eq.0) then

c           Unsymmetric tangent: Ibrahimbegovic Incremental Rotation

            if(nint(d(17)).eq.2 .or. nint(d(17)).eq.3) then
              call framf3d(d,ul,xl,ix,s,r,ndf,ndm,nst,isw)

c           Energy conserving:

            elseif(nint(d(17)).eq.4) then
              call framf3e(d,ul,xl,ix,s,r,ndf,ndm,nst,isw)

c           Unsymmetric tangent: Co-rotational

            elseif(nint(d(17)).eq.6) then
              call framf3c(d,ul,xl,s,r,ndf,ndm,nst,isw)

c           Unsymmetric tangent: Simo-VuQuoc Displacement

            else
              call framf3b(d,ul,xl,s,r,ndf,ndm,nst,isw)

            endif

c         No shear: 2-node, cubic, second order strains

          else
            call franf3d(d,ul,xl,s,r,ndf,ndm,nst,isw)
          endif

        endif

        if(isw.eq.3) then
          ctan(1) = ctan1
        endif

      endif

c     Format statements

2000  format(5x,'T h r e e   D i m e n s i o n a l   F r a m e'/)

2001  format(/5x,'*ERROR* Some of following values are wrong:'
     &       /6x,'ndf(should be 6) =',i3/6x,'ndm (should be 3) =',i3)

2002  format(/5x,'*ERROR* No reference vector definition for'
     &          ,' axes: lref =',i3)

2003  format(/5x,'*ERROR* Cannot use CROSs SECTion data'
     &          ,' with an inelastic material model.'/
     &       13x,'Inelastic material model requires SECTion type.')
      end
