c$Id:$
      subroutine presld(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set pstyp = -1 for no plots                      02/01/2007
c          Add capability for 8 and 9-node elements
c       2. Replace d(182) by inputs of 'noda' and 'gaus'    20/01/2007
c       3. Set pstyp =  0 for no plots                      31/08/2008
c       4. Correct axisymmetric computation for follower
c          load computations to account for radius change.  20/04/2011
c          Add output of geometry type.  Allow for plot of
c          pressure mesh ('plot','nopl' option).
c       5. Compute reactions for pressure load (remove the  17/11/2011
c          .not.ddfl). Remove 'c_tanfl.h' include.
c       6. Set pt1 and pt2 to 0.0 if plane problem          07/12/2012
c       7. Change shape function from 1-d to 2-d for height 10/12/2012
c          computation of 3-d problems. Change sg2 to sg1,
c          sg3 to sg2, shp2 to shp1 and shp3 to shp2.
c       8. Use 'quadr2d' and 'interp2d' to get triangular   11/12/2012
c          facets also. Set d(182), d(5) for quadrature.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Nodal force and tangent array for pressure loading
c                Includes case for dead load as well as follower
c                pressure loading

c      INPUT variables
c        d(60)      Value of constant pressure on face
c        xl(4,*)    Nodal coordinates
c        ul(4,*)    Nodal displacements
c        ndf        Number of DOF / node
c        ndm        Space dimension
c        nst        Dimension of residual vector

c      OUTPUT variables
c        p(ndf,4)   Contribution to residual
c        s(nst,nst) Contribution to striffness matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'eltran.h'
      include  'iofile.h'
      include  'pglob1.h'
      include  'prld1.h'
      include  'qudshp.h'

      logical   pcomp, errck, tinput, mread, prsfl, nodfl
      character texti*15
      integer   l, ndf, ndm, nst, isw, ii, jj, kk, i1, j1, ix(*)
      real*8    pn, pp, pr, pt1,pt2, x0,v0, xx, rr
      real*8    td(7), xu(3,9),dx(3,2)
      real*8    d(*), xl(ndm,*), ul(ndf,*), p(ndf,*), s(nst,*)

      save

      if(isw.eq.1) then

c       Global axisymmetric type indicator

        if(ndm.eq.2) then
          d(81) = g2type
        endif

c       Global check for finite deformation (follower loading)

        if(gdeffl) then
          d(82) = 1.0d0
        else
          d(82) = 0.0d0
        endif

c       Deactivate dof in element for dof > ndm

        do ii = ndm+1,ndf
          ix(ii) = 0
        end do ! ii

c       Input data for pressure loading

        nodfl = .false.
        mread = .true.
        do while (mread)
          errck = tinput(texti,1,td,7)
          if(pcomp(texti,'load',4) .or. pcomp(texti,'pres',4)) then
            prsfl = .true.
            d(60) = td(1)
            d(61) = nint(td(2))
          elseif(pcomp(texti,'grad',4) .or. pcomp(texti,'bend',4)) then
            d(62) = td(1)
            d(63) = nint(td(2))
          elseif(pcomp(texti,'hydr',4)) then
            prsfl = .false.
            d(64) = td(1) ! j   = Direction of pressure gradient
            d(65) = td(2) ! x_j = Zero pressure elevation
            d(66) = td(3) ! v_j = Elevation velocity
            d(67) = td(4) ! rho = Density
            d(68) = td(5) ! g_j = Gravity value
            d(69) = td(6) ! p_n = Proportional load number pressure
            d(70) = td(7) ! p_h = Proportional load number height
          elseif(pcomp(texti,'axis',4)) then
            d(81) = 3.0d0
          elseif(pcomp(texti,'plan',4)) then
            d(81) = 0.0d0
          elseif(pcomp(texti,'fini',4) .or. pcomp(texti,'foll',4)) then
            d(82) = 1.d0
          elseif(pcomp(texti,'smal',4) .or. pcomp(texti,'dead',4)) then
            d(82) = 0.d0
          elseif(pcomp(texti,'plot',4)) then   ! Allows plot of pressure
            d(85) = 1.d0                       ! Mesh
          elseif(pcomp(texti,'nopl',4)) then   ! Turns off mesh plot
            d(85) = 0.d0
          elseif(pcomp(texti,'quad',4)) then   ! Nodal quadrature
            d( 5) = min(5.d0,td(1))
          elseif(pcomp(texti,'noda',4)) then   ! Nodal quadrature
            d(182) = 1.d0
            nodfl =.true.
          elseif(pcomp(texti,'gaus',4)) then   ! Gauss-Legendre quadr
            d(182) = 1.d0
            nodfl =.false.
          else
            mread = .false.
          endif
        end do ! while

c       Output input definitions

        if(prsfl) then     ! Pressure: Constant and gradient
          write(iow,2000) d(60),nint(d(61)),d(62),nint(d(63))
          if(d(82).gt.0.0d0) then
            write(iow,2002)
          endif
          d(83) = 0.0d0
        else
          write(iow,2001) nint(d(64)),d(65),d(66),d(67),d(68),
     &                    nint(d(69)),nint(d(70))
          d(83) = 1.0d0
        endif
        if(ndm.eq.2) then
          if(nint(d(81)).ne.3) then
            write(iow,2005) 'Plane'
          elseif(nint(d(81)).eq.3) then
            write(iow,2005) 'Axisymmetric'
          endif
        elseif(ndm.eq.3) then
          write(iow,2005) 'Three dimensional'
        endif

c       Output quadrature type

        if(nodfl) then
          write(iow,2003)
        else
          write(iow,2004)
        endif

c       Set for plots if requested (For testing mesh validity)

        if(nint(d(85)).eq.0) then
          pstyp = 0
        else
          pstyp = ndm - 1
        endif

c     Compute residual and tangent

      elseif(isw.eq.3 .or. isw.eq.6 .or. isw.eq.23) then

c       Compute nodal coordinates in correct reference frame

        if(d(82).eq.0.0d0) then ! Fixed force  (small option)
          do ii = 1,nel
            do jj = 1,ndm
              xu(jj,ii) = xl(jj,ii)
            end do ! jj
          end do ! ii
        else                    ! Follower force  (finite option)
          do ii = 1,nel
            do jj = 1,ndm
              xu(jj,ii) = xl(jj,ii) + ul(jj,ii)
            end do ! ii
          end do ! ii
        endif

c       2-D Problems

        if(ndm.eq.2) then

c         Get quadrature information

          lint = nel
          if(nint(d(182)).gt.0) then  ! Nodal quadrature
            call int1dn(lint, sg1)
          else                       ! Gauss quadrature
            call int1d (lint, sg1)
          endif

c         Loop over quadrature points

          do l = 1,lint

c           Compute shape functions & parent derivatives

            call shap1d( sg1(1,l), nel, shp1(1,1,l) )

c           Compute coordinate gradient w/r parent coordinates

            do ii = 1,ndm
              dx(ii,1) = shp1(1,1,l)*xu(ii,1)
              do jj = 2,nel
                dx(ii,1) = dx(ii,1) + shp1(1,jj,l)*xu(ii,jj)
              end do ! jj
            end do ! ii

c           Compute radius for axisymmetric problem

            if(nint(d(81)).eq.3) then
              rr = 0.0d0
              do jj = 1,nel
                rr = rr + shp1(2,jj,l)*xu(1,jj)
              end do ! jj
            else
              rr = 1.0d0         ! Set radius = 1 for plane/3-d problems
            endif

c           Compute nodal loads for pressures

            if(nint(d(83)).eq.0) then          ! Uniform pressure
              kk = nint(d(61))
              if(kk.eq.0) then
                pn = dm*d(60)*sg1(2,l)         ! dm = total prop ld.
              else
                pn = prldv(kk)*d(60)*sg1(2,l)  ! For specified prop ld.
              endif
              kk = nint(d(63))                 ! Add gradient term
              if(kk.eq.0) then
                pn = pn + dm*d(62)*sg1(1,l)*sg1(2,l)
              else
                pn = pn + prldv(kk)*d(62)*sg1(1,l)*sg1(2,l)
              endif
            elseif(nint(d(83)).eq.1) then      ! ??
              ii = nint(d(64))
              x0 = d(65)
              v0 = d(66)
              kk = nint(d(70))
              if(kk.eq.0) then
                x0 = x0 + v0*dm
              else
                x0 = x0 + v0*prldv(kk)
              endif
              xx = 0.0
              do jj = 1,nel
                xx = xx + shp1(2,jj,l)*xl(ii,jj)
              end do ! jj
              if(xx.le.x0) then
                pn = d(67)*d(68)*(xx - x0)*sg1(2,l)
                kk = nint(d(69))
                if(kk.eq.0) then
                  pn = pn*dm
                else
                  pn = pn*prldv(kk)
                endif
              else
                pn = 0.0d0
              endif
            endif

c           Set surface pressure x radius

            pr = pn*rr                   ! Axisymmetric correction

c           Compute residual for pressure

            do ii = 1,nel
              pp      = shp1(2,ii,l)*pr
              p(1,ii) = p(1,ii) + pp*dx(2,1)
              p(2,ii) = p(2,ii) - pp*dx(1,1)
            end do ! ii

c           Compute tangent for follower load case

            if(d(82).gt.0.0d0) then
              i1 = 0
              do ii = 1,nel
                pp  = pn*shp1(2,ii,l)*ctan(1)
                if(nint(d(81)).eq.3) then
                  pt1 = pp*dx(2,1)
                  pt2 = pp*dx(1,1)
                else
                  pt1 = 0.0d0
                  pt2 = 0.0d0
                endif
                pp  = pp*rr
                j1  = 0
                do jj = 1,nel
                  s(i1+1,j1+1) = s(i1+1,j1+1) - pt1*shp1(2,jj,l)
                  s(i1+1,j1+2) = s(i1+1,j1+2) - pp *shp1(1,jj,l)

                  s(i1+2,j1+1) = s(i1+2,j1+1) + pt2*shp1(2,jj,l)
     &                                        + pp *shp1(1,jj,l)
                  j1 = j1 + ndf
                end do ! jj
                i1 = i1 + ndf
              end do ! ii
            endif

          end do ! l

c       3-D Problems

        elseif(ndm.eq.3) then

c         Get quadrature information

          call quadr2d(d,.false.) !Returns quadrature for various shapes

c         Loop over quadrature points

          do l = 1,lint

c           Compute shape functions and geometric factors
c           The .true. gives derivatives w/r parent coordinates

            call interp2d(l,xl,ix,ndm,nel,.true.)

c           Compute tangent vectors in 2 parent coordinate distances

            do ii = 1,3
              dx(ii,1) = 0.0d0
              dx(ii,2) = 0.0d0
              do jj = 1,nel
                dx(ii,1) = dx(ii,1) + shp2(1,jj,l)*xu(ii,jj)
                dx(ii,2) = dx(ii,2) + shp2(2,jj,l)*xu(ii,jj)
              end do ! jj
            end do ! ii

c           Compute nodal loads for pressures

            if(nint(d(83)).eq.0) then
              kk = nint(d(61))
              pn = d(60)*sg2(3,l)
              if(kk.eq.0) then
                pn = pn*dm
              elseif(kk.gt.0) then
                pn = pn*prldv(kk)
              endif
            elseif(nint(d(83)).eq.1) then
              ii = nint(d(64))
              x0 = d(65)
              v0 = d(66)
              kk = nint(d(70))
              if(kk.eq.0) then
                x0 = x0 + v0*dm
              else
                x0 = x0 + v0*prldv(kk)
              endif
              xx = 0.0
              do jj = 1,nel
                xx = xx + shp2(3,jj,l)*xl(ii,jj)
              end do ! jj
              if(xx.le.x0) then
                pn = d(67)*d(68)*(xx - x0)*sg2(3,l)
                kk = nint(d(69))
                if(kk.eq.0) then
                  pn = pn*dm
                else
                  pn = pn*prldv(kk)
                endif
              else
                pn = 0.0d0
              endif
            endif

c           Compute residual

            do ii = 1,nel
              pp      = shp2(3,ii,l)*pn
              p(1,ii) = p(1,ii) + pp*(dx(2,1)*dx(3,2) - dx(3,1)*dx(2,2))
              p(2,ii) = p(2,ii) + pp*(dx(3,1)*dx(1,2) - dx(1,1)*dx(3,2))
              p(3,ii) = p(3,ii) + pp*(dx(1,1)*dx(2,2) - dx(2,1)*dx(1,2))
            end do ! ii

c           Compute tangent for follower load

            if(d(82).gt.0.0d0) then
              i1 = 0
              do ii = 1,nel
                pp = shp2(3,ii,l)*pn*ctan(1)
                j1 = 0
                do jj = 1,nel
                 s(i1+1,j1+2) = s(i1+1,j1+2) - pp*(shp2(1,jj,l)*dx(3,2)
     &                                       -     shp2(2,jj,l)*dx(3,1))
                 s(i1+2,j1+1) = s(i1+2,j1+1) + pp*(shp2(1,jj,l)*dx(3,2)
     &                                       -     shp2(2,jj,l)*dx(3,1))

                 s(i1+2,j1+3) = s(i1+2,j1+3) - pp*(shp2(1,jj,l)*dx(1,2)
     &                                       -     shp2(2,jj,l)*dx(1,1))
                 s(i1+3,j1+2) = s(i1+3,j1+2) + pp*(shp2(1,jj,l)*dx(1,2)
     &                                       -     shp2(2,jj,l)*dx(1,1))

                 s(i1+3,j1+1) = s(i1+3,j1+1) - pp*(shp2(1,jj,l)*dx(2,2)
     &                                       -     shp2(2,jj,l)*dx(2,1))
                 s(i1+1,j1+3) = s(i1+1,j1+3) + pp*(shp2(1,jj,l)*dx(2,2)
     &                                       -     shp2(2,jj,l)*dx(2,1))
                 j1 = j1 + ndf
                end do ! jj
                i1 = i1 + ndf
              end do ! ii
            endif

          end do ! l

        endif ! ndm test

c     External node check

      elseif(isw.eq.26) then

      endif ! isw test

c     Formats

2000  format(5x,'P r e s s u r e   L o a d i n g'//
     &      10x,'Loading Intensity    ',1p,1e12.4/
     &      10x,'Proportional Load No.',i12/
     &      10x,'Gradient Intensity   ',1p,1e12.4/
     &      10x,'Proportional Load No.',i12)

2001  format(5x,'P r e s s u r e   L o a d i n g'//
     &      10x,'Loading direction         =',i12/
     &      10x,'Zero pressure elevation   =',1p,1e12.4/
     &      10x,'Velocity of surface       =',1p,1e12.4/
     &      10x,'Density/specific weight   =',1p,1e12.4/
     &      10x,'Gravity acceleration      =',1p,1e12.4/
     &      10x,'Proportional Load Pressure=',i12/
     &      10x,'Proportional Load Height  =',i12)

2002  format(10x,'Follower Loading')
2003  format(10x,'Nodal quadrature')
2004  format(10x,'Gauss quadrature')
2005  format(10x,a,' analysis')

      end
