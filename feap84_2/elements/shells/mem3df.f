c$Id:$
      subroutine mem3df ( d, ul, xl, s, p, ndf, ndm, nst, isw )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Description:    MEM3DF is a 4-node bilinear element for
c                     static/transient finite deformation analysis
c                     of 3-dimensional elastic membrane problems.

c     Date:           July 1998

c     Author:         R.L. Taylor (Based on shell element by
c                                  J.C. Simo, et. al.)

c     Control Information:
c     -------------------------
c     ndm = 3 ....... X, Y, Z Cartesian Coordinates at Nodes.
c     ndf = 3 ....... U-x, U-y, U-z, DOF at Nodes.
c     nen =>3 ....... # of Element Nodes.
c     nst = 12 ...... # of Element DOF's.

c     Element Material Set Input:
c     --------------
c     mate, Material-Group-Number
c       membrane
c       Elastic Isotropic  E, Nu (or Fung Exponential Modal)
c       Thickness membrane thk
c       Density   mass     Rho
c       finite

c     Input Variables:
c     ----------------
c     E ............. Young's modulus.
c     Nu ............ Poisson's ratio.
c     thk ........... Membrane thickness. If set to zero, membrane
c                     thickness is assumed to be average of
c                     global thicknesses defined at nodes.
c                     If non-zero, this value overrides global
c                     input.
c     Rho ........... Translational density.

c     FEAP Common Arrays:
c     -------------------
c     hr............. History variable block.

c     FEAP Solution Array:
c     --------------------
c     ul ............ Localized displacements and Velocities:
c                     ul(5,Node,1) := Total displacement T0-Tn+a.
c                     ul(5,Node,2) := Step displacement Tn-Tn+a.
c                     ul(5,Node,3) := Displacement increment. Tn+1
c                     ul(5,Node,4) := Velocity at Tn+a.
c                     ul(5,Node,5) := Acceleration at Tn+a.

c     Major Global Variables:
c     -----------------------
c     thkl .......... Localized nodal thickness.

c     Major Element Variables:
c     ------------------------
c     shp ........... Nodal shape function values.
c     shp1,shp2 ..... Nodal shape function natural derivatives.
c     shx1,shx2 ..... Nodal shape function cartesian derivatives.
c     sg ............ Gauss point coordinates and weights.
c     b ............. Discrete strain-displacement operator.
c     c ............. Resultant stress constitutive matrix.
c     zphG .......... Initial local coordinates
c     zph1,zph2 ..... Initial coordinate local derivatives
c     zpx1,zpx2 ..... Initial coordinate global derivatives
c     cphi .......... Current local nodal coordinates.
c     cpx1,cpx2 ..... Current coordinate global derivatives
c     xjinv ......... Local-global coordinate Jacobian matrix.
c     xjs ........... Reference jacobian at Gauss points.
c     xjw ........... Reference jacobian weighted for integration.
c     ze,ce ......... Initial, current membrane strain measure.
c     sn ............ Membrane stress.
c     aa ............ Gauss point accelerations at t-n+a.
c     nu ............ Poisson's ratio.
c     rhoa .......... Trans. thickness weighted density multiplyied
c                     by jacobian and weight at Gauss point.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'ddata.h'
      include  'crotas.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'erotas.h'
      include  'fdata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'prstrs.h'
      include  'strnum.h'
      include  'tdata.h'
      include  'comblk.h'

      logical   rhs , lhs , errck
      integer   i,j,k,l, lint, i5,i3,j5,j3, ndf,ndm,nst,isw
      real*8    rhoa , xjdet , xlgnrm , xj1,xj2,xj3, dthk,fac1,facn
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),s(nst,*),p(nst)
      real*8    shp(4,4),shp1(4,4),shp2(4,4),shx1(4,4),shx2(4,4)
      real*8    sg(3,4),fphi(3,4),bphi(3,4),zph1(3,4),zph2(3,4)
      real*8    zpx1(3,4),zpx2(3,4),cpx1(3,4),cpx2(3,4),cphi(3,4)
      real*8    bpx1(3,4),bpx2(3,4),fpx1(3,4),fpx2(3,4),ze(3,4)
      real*8    be(3,4),ce(3,4),fe(3,4),sn(3,4),aa(3)
      real*8    ploc(12),sloc(12,12),c(3,3),b(3,12,4),b1(3,12,4)
      real*8    xjinv(2,2),xjs(4),xjw(4),xlg (3,3)

      save

c     Set Thickness From Global Input

      if (isw.gt.1 .and. d(14).lt.1.d-5) then
        dthk  = 0.25d0 * (thkl(1)+thkl(2)+thkl(3)+thkl(4))
        if(dthk.lt.0.0d0) then
          write(iow,3002) dthk
          stop ' Membrane element error'
        endif
      else
        dthk  = d(14)
      endif

c     Input Material Properties:

      if(isw.eq.1) then

c       Check input data

        errck = .false.
        if(nint(d(20)).eq.7) then
          if(d(30).le.0.0d0) then
            write(iow,3000)
            errck = .true.
          endif
        else
          if(d(1).le.0.0d0) then
            write(iow,3000)
            errck = .true.
          endif
        endif
        if(d(14).le.0.0d0) then
          write(iow,3001)
        endif
        if(errck) stop '  Membrane material property errors'

c       Initialize Shape Functions.

        l     = 2
        call int2d      ( l        , lint    , sg )
        do i = 1 , lint
          call me3fshap ( sg (1,i), shp(1,i), shp1(1,i), shp2(1,i)  )
        end do ! i

c     Compute Element Dynamic Tangent and Residual Array:

      elseif(isw.ge.3 .and. isw.le.6  .or. isw.eq.8 .or. isw.eq.13) then

c       Set points

c       Compute Midsurface Gauss-Point Interpolation

        do i = 1 , lint
          do j = 1 , 3

c           Midsurface

            zph1(j,i) = shp1(1,i)*xl(j,1) + shp1(2,i)*xl(j,2)
     &                + shp1(3,i)*xl(j,3) + shp1(4,i)*xl(j,4)
            zph2(j,i) = shp2(1,i)*xl(j,1) + shp2(2,i)*xl(j,2)
     &                + shp2(3,i)*xl(j,3) + shp2(4,i)*xl(j,4)
          end do ! j

c         Compute Reference Surface Jacobian

          xj1    =  zph1(2,i)*zph2(3,i) - zph1(3,i)*zph2(2,i)
          xj2    = -zph1(1,i)*zph2(3,i) + zph1(3,i)*zph2(1,i)
          xj3    =  zph1(1,i)*zph2(2,i) - zph1(2,i)*zph2(1,i)
          xjs(i) =  sqrt( xj1**2 + xj2**2 + xj3**2 )
          xjw(i) =  sg(3,i) * xjs(i)

        end do ! i

c       Do mass computation

        if(isw.eq.5) then
          call me3fmasg( d, dthk, shp, xjw, s, lint, ndf, nst)
          return
        endif

c       Set Flags From Input

        lhs    = isw.eq.3
        rhs    = .true.

c       Compute Current Configuration

        fac1 = 1.d0 / theta(3) - 1.d0
        do i = 1 , 4
          do j = 1 , 3
            cphi(j,i) = xl(j,i)   + ul(j,i,1)
            bphi(j,i) = cphi(j,i) - ul(j,i,2)
            fphi(j,i) = cphi(j,i) + ul(j,i,2)*fac1
          end do ! j
        end do ! i

c       Compute Gauss-Point Interpolations and Derivatives

        fac1 = theta(3)
        facn = 1.d0 - fac1

        do i = 1 , lint

c         Compute Surface Normal

          call vecp ( zph1(1,i) , zph2(1,i) , xlg(1,3) )

          xlgnrm   = 1.d0/sqrt(xlg(1,3)**2 + xlg(2,3)**2 + xlg(3,3)**2)
          xlg(1,3) = xlg(1,3) * xlgnrm
          xlg(2,3) = xlg(2,3) * xlgnrm
          xlg(3,3) = xlg(3,3) * xlgnrm

c         Compute Local-Global Jacobian

c         call me3flmda ( xlg(1,3) , xlg  )
          call me3flmda ( zph1(1,i) , zph2(1,i) , xlg  )

          xjinv(2,2) = xlg(1,1)*zph1(1,i)
     &               + xlg(2,1)*zph1(2,i)
     &               + xlg(3,1)*zph1(3,i)
          xjinv(2,1) = xlg(1,2)*zph1(1,i)
     &               + xlg(2,2)*zph1(2,i)
     &               + xlg(3,2)*zph1(3,i)
          xjinv(1,2) = xlg(1,1)*zph2(1,i)
     &               + xlg(2,1)*zph2(2,i)
     &               + xlg(3,1)*zph2(3,i)
          xjinv(1,1) = xlg(1,2)*zph2(1,i)
     &               + xlg(2,2)*zph2(2,i)
     &               + xlg(3,2)*zph2(3,i)

          xjdet      = 1.d0/(xjinv(1,1)*xjinv(2,2)
     &                     - xjinv(1,2)*xjinv(2,1))

          xjinv(1,1) = xjinv(1,1) * xjdet
          xjinv(1,2) =-xjinv(1,2) * xjdet
          xjinv(2,1) =-xjinv(2,1) * xjdet
          xjinv(2,2) = xjinv(2,2) * xjdet

c         Compute Global Shape Function Derivatives

          call me3fmdsh ( shp1(1,i) , shp2(1,i) , xjinv     ,
     &                    shx1(1,i) , shx2(1,i) )

c         Interpolate Global Derivatives

          do j = 1 , 3

c           Midsurface

            zpx1(j,i) = shx1(1,i)*xl  (j,1)    + shx1(2,i)*xl  (j,2)
     &                + shx1(3,i)*xl  (j,3)    + shx1(4,i)*xl  (j,4)
            zpx2(j,i) = shx2(1,i)*xl  (j,1)    + shx2(2,i)*xl  (j,2)
     &                + shx2(3,i)*xl  (j,3)    + shx2(4,i)*xl  (j,4)
            cpx1(j,i) = shx1(1,i)*cphi(j,1)    + shx1(2,i)*cphi(j,2)
     &                + shx1(3,i)*cphi(j,3)    + shx1(4,i)*cphi(j,4)
            cpx2(j,i) = shx2(1,i)*cphi(j,1)    + shx2(2,i)*cphi(j,2)
     &                + shx2(3,i)*cphi(j,3)    + shx2(4,i)*cphi(j,4)
            fpx1(j,i) = shx1(1,i)*fphi(j,1)    + shx1(2,i)*fphi(j,2)
     &                + shx1(3,i)*fphi(j,3)    + shx1(4,i)*fphi(j,4)
            fpx2(j,i) = shx2(1,i)*fphi(j,1)    + shx2(2,i)*fphi(j,2)
     &                + shx2(3,i)*fphi(j,3)    + shx2(4,i)*fphi(j,4)
            bpx1(j,i) = shx1(1,i)*bphi(j,1)    + shx1(2,i)*bphi(j,2)
     &                + shx1(3,i)*bphi(j,3)    + shx1(4,i)*bphi(j,4)
            bpx2(j,i) = shx2(1,i)*bphi(j,1)    + shx2(2,i)*bphi(j,2)
     &                + shx2(3,i)*bphi(j,3)    + shx2(4,i)*bphi(j,4)

          end do ! j

c         Compute Strain Measures

          call me3fstrn ( zpx1(1,i)     , zpx2(1,i)   , ze  (1,i)  )
          call me3fstrn ( fpx1(1,i)     , fpx2(1,i)   , fe  (1,i)  )

          if (isw.eq.13 .or. .not.fl(9)) then
            do j = 1,3
              ce(j,i) = fe(j,i) - ze(j,i)
            end do ! j
          else
            call me3fstrn ( bpx1(1,i)     , bpx2(1,i)   , be  (1,i)  )
            do j = 1 , 3
              ce(j,i) = fac1*fe(j,i) + facn*be(j,i) - ze(j,i)
            end do ! j
          endif

c         Set Up Strain-Displacement Matrices

          call me3fbmtx ( shx1(1,i)     , shx2(1,i)   ,
     &                    cpx1(1,i)     , cpx2(1,i)   ,
     &                    b(1,1,i)                    )

          if (lhs) then

            call me3fbmtx ( shx1(1,i)     , shx2(1,i)   ,
     &                      fpx1(1,i)     , fpx2(1,i)   ,
     &                      b1(1,1,i)                   )
          endif

        end do ! i

c       Zero full (local) stiffness and residuals

        do i = 1,12
          do j = 1,12
            sloc(j,i) = 0.0d0
          end do ! j
          ploc(i) = 0.0d0
        end do ! i

c       Galerkin and Shear Treatment

        do i = 1,lint

c         Compute Constitutive Relations (weighted by jacobian)

          call me3fcmtx ( d , ce(1,i) , dthk  , xjw(i) , c )

c         Compute Stresses (weighted by jacobian)

          call me3fstrs ( d , c  , ce(1,i) , dthk, xjw(i), sn(1,i)  )

c         Store time history plot data for element

          i3 = 3*(l-1)
          do j = 1,3
            tt(j+i3) = sn(j,i)
          end do ! j

c         Compute Residual

          call me3fresd ( b(1,1,i) , sn(1,i) , ploc )

          if (lhs) then

c           Compute Material Stiffness

            call me3fstif ( b(1,1,i) , b1(1,1,i) , c , sloc )

c           Compute Geometric Stiffness

            call me3fgeom ( shx1(1,i) , shx2(1,i) , sn(1,i) , c , sloc )
          endif

        end do ! i

c       Compute Transient Terms

        if (fl(9) .and. noi.ne.6) then

          do i = 1,lint

c           Interpolate Accelerations

            do j = 1 , 3
              aa(j) = shp(1,i)*ul(j,1,5) + shp(2,i)*ul(j,2,5)
     &              + shp(3,i)*ul(j,3,5) + shp(4,i)*ul(j,4,5)
            end do ! j

c           Compute Transient Residual and Tangent

            rhoa = d(4) * d(14) * xjw(i)
            call me3ftran ( dt,    shp(1,i)  , aa     ,
     &                      rhoa,  ploc      , sloc )

          end do ! l

        endif

c       Stress and energy processes

        if (isw.eq.4 .or. isw.eq.8)  then
          call me3fstre ( xl, ce, sn, shp, xjs, lint, ndm, isw )
          rhs = .false.
        elseif (isw.eq.13) then
          call me3fener ( d, dthk, ul, shp, fphi, xjw, ndf, lint )
          rhs = .false.
        elseif(lhs) then

c         Scale tangent stiffness for total increment du_t-n+1

          do j = 1,12
            do i = 1,12
              sloc(i,j) = fac1*sloc(i,j)
            end do ! i
          end do ! j

c         Transform (Reduce) Stiffness/Residual - First Node Loop

          i3 = 0
          i5 = 0
          do i = 1 , 4

            j3 = 0
            j5 = 0
            do j = 1 , 4

              do l = 1 , 3
                do k = 1 , 3
                  s(i5+k,j5+l) = sloc(i3+k,j3+l)
                end do ! k
              end do ! l

              j3 = j3 + 3
              j5 = j5 + ndf
            end do ! j

            i3 = i3 + 3
            i5 = i5 + ndf
          end do ! i

        endif

        if(rhs) then

          i3 = 0
          i5 = 0
          do i = 1 , 4

            do j = 1 , 3
              p(i5+j) = ploc(i3+j)
            end do ! j

            i3 = i3 + 3
            i5 = i5 + ndf
          end do ! i

c         Modify residual and stiffness for pressure & body loading

          call me3fsurf(d,xl,ul,xjw, ndf,ndm,nst, p,s)

        endif
      endif

c     Formats:

3000  format(' *ERROR* Modulus is zero or negative for membrane')

3001  format(' *WARNING* Thickness zero: Compute from nodal input.')

3002  format(' *ERROR* Thickness is zero or negative for membrane')

      end
