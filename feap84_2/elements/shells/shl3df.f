c$Id:$
      subroutine shl3df ( d, ul, xl, s, p, ndf, ndm, nst, isw )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add check on noi > 0 for transient analysis      01/05/2012
c       2. Add dthk on call to sh3fsurf                     31/05/2013
c       3. Set plot weight xjw(i) to 1 for projections      01/06/2013
c-----[--.----+----.----+----.-----------------------------------------]
c       Description:    SHL3DF is a 4-node bilinear element for
c                       static/transient finite deformation analysis
c                       of inextensible 3-dimensional elastic shell
c                       problems.

c       Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c       Date:           January 1991.

c       Revised:        R.L. Taylor  - - February 1997

c       Modification:   J.C. Simo and N. Tarnow to handle
c                       Energy-Momentum method and shell-intersections.
c                       January 1992.
c                       R.L. Taylor, Adapted to version 6.0 and later.
c                       february 1997.

c       features:       Static and transient inextensible analysis.
c                       Full Isoparametric interpolation, including
c                             inextensible director field.
c                       2-DOF director increments (Delta-T).
c                       Bathe-Dvorkin shear treatment.
c                       Displacement version or
c                       Pian-Summihara membrane/bending treatment.

c       Caution:        This element cannot be used with rotational
c                       increments in excess of 180-degrees.

c       References:     J.C. Simo, D.D. Fox, M.S. Rifai,
c                       "On a Stress Resultant Geometrically Exact
c                       Shell Model. Parts I-VII,"
c                       Comp. Meth. Appl. Mech. & Engng, 1989-1991.
c-----[--.----+----.----+----.-----------------------------------------]
c       Control Information:
c       -------------------------
c       ndm = 3 ....... X, Y, Z Cartesian Coordinates at Nodes.
c       ndf = 5 ....... U-x, U-y, U-z, DT-1, DT-2 DOF at Nodes.
c       nen = 4 ....... # of Element Nodes: Input Counterclockwise.
c       nst = 20 ...... # of Element DOF's.
c-----[--.----+----.----+----.-----------------------------------------]
c       Element Material Set Input:
c       --------------
c       mate, Material-Group-Number
c         shell
c         Elastic Isotropic  E, Nu
c         Thickness shell    h, kappa
c         Density   mass     Rho, rot-mass-factor
c         displacement or mixed
c         finite

c       Input Variables:
c       ----------------
c       E ............. Young's modulus.
c       Nu ............ Poisson's ratio.
c       G-s ........... Transverse shear modulus - including
c                       shear correction factor.
c       thk ........... Shell thickness. If set to zero, shell
c                       thickness is assumed to be average of
c                       global thicknesses defined at nodes.
c                       If non-zero, this value overrides global
c                       input.
c       Rho-trn ....... Translational density.
c       Rho-rot ....... Rotational density.
c                       The reason for inputing two different density
c                       values is to be able to run fictitious
c                       material properties invented in
c                       Simo & Vu-Quoc beam formulation papers,
c                       e.g., pencil toss problem.
c       Opt-mix ....... Mixed treatment option:
c                       0 => Galerkin (displacement formulation),
c                       1 => Hellinger-Reissner (Pian-Summihara).
c-----[--.----+----.----+----.-----------------------------------------]

c       Quasi-Static Element Operation:
c       -------------------------------
c       DT,, Time-Step-Size
c       TOL,, Tolerance-Level           (Dflt: 1.d-10, Sugstd: 1.d-15)
c       PROP,, Proportional-Load-Option (Dflt: 1 - Continuous fctn)
c           Proportional-Load-Input     (Dflt: Linear with time)
c       LOOP,, Number-Of-Time-Steps
c         TIME
c         LOOP,, Number-Of-Newton-Iterations
c           TANG,,1
c         NEXT
c       NEXT

c       Transient Element Operation:
c       ----------------------------
c       DT,, Time-Step-Size
c       TOL,, Tolerance-Level           (Dflt: 1.d-10, Sugstd: 1.d-15)
c       PROP,, Proportional-Load-Option (Dflt: 1 - Continuous fctn)
c           Proportional-Load-Input     (Dflt: Linear with time)
c       TRAN, CONS, Beta, Gamma, Alpha  (Dflt: 0.5, 1.0, 0.5)
c       CMAS
c       LOOP,, Number-Of-Time-Steps
c         TIME
c         LOOP,, Number-Of-Newton-Iterations
c           TANG,,1
c         NEXT
c       NEXT

c       ArcLength Operation:
c       --------------------
c       DT,, 1.0                        (Time step must be chosen 1)
c       TOL,, Tolerance-Level           (Dflt: 1.d-10, Sugstd: 1.d-15)
c       PROP,, Proportional-Load-Option (Dflt: 1 - Continuous fctn)
c           Proportional-Load-Input     (Dflt: Linear with time)
c       ARCL,, 2
c       LOOP,, Number-Of-Time-Steps
c           TIME
c           LOOP,, Number-Of-Newton-Iterations
c               TANG,,1
c           NEXT
c       NEXT

c       Note:
c       -----
c       ARCL,, 1 ...... Can be used to switch arclength direction
c                       during execution.

c-----[--.----+----.----+----.-----------------------------------------]
c       FEAP Common Arrays:
c       -------------------
c       hr............. History variable block.
c       mh ............ Global FEAP block common.

c       FEAP Solution Array:
c       --------------------
c       ul ............ Localized displacements and Velocities:
c                       ul(5,Node,1) := Total displacement T0-Tn+a.
c                       ul(5,Node,2) := Step displacement Tn-Tn+a.
c                       ul(5,Node,3) := Displacement increment. Tn+1
c                       ul(5,Node,4) := Velocity at Tn+a.
c                       ul(5,Node,5) := Acceleration at Tn+a.

c       Major Global Variables:
c       -----------------------
c       enflag ........ Energy conservation flag.
c       etotn,etot1 ... Total Energy at t-n and t-n+1.
c       epotn,epot1 ... Potential Energy at t-n and t-n+1.
c       alphar ........ Alpha residual.
c       alphat ........ Alpha tangent.
c       alphak ........ Alpha at previous iteration.
c       dalpha ........ Alpha Increment.
c       xln ........... Localized nodal orthogonal transformation
c                       matrices, xln(3,3,ElementNode,Ntime):
c                       Ntime=1 => Value at time t-n,
c                       Ntime=2 => Value at time t-n+alpha,
c                       Ntime=3 => Value at time t-n+1,
c                       Ntime=4 => Value at time 0.
c       rots .......... Localized nodal rotation vectors,
c                       rots(3,ElementNode,1) = Incremental rotation,
c                       rots(3,ElementNode,2) = Total rotation t-n+1.
c       rvel .......... Localized nodal rot. velocity vectors,
c                       rvel(3,ElementNode,Ntime):
c                       Ntime=1 => Value at time t-n,
c                       Ntime=2 => Value at time t-n+1.
c       racc .......... Localized nodal rot. acceleration vectors,
c                       racc(3,ElementNode,Ntime):
c                       Ntime=1 => Value at time t-n,
c                       Ntime=2 => Value at time t-n+1.
c       thkl .......... Localized nodal thickness.
c-----[--.----+----.----+----.-----------------------------------------]
c       Major Element Variables:
c       ------------------------
c       m ............. Transformation matrix current index:
c                       In stiffness assembly, m=2 => t-n+alpha,
c                       In time-step output,   m=3 => t-n+1.
c       dir ........... Element local directors
c                       dir(3,ElementNode,Ntime):
c                       Ntime=1 => Value at time t-n,
c                       Ntime=2 => Value at time t-n+alpha,
c                       Ntime=3 => Value at time t-n+1,
c                       Ntime=4 => Value at time t-0,
c                       Ntime=5 => Value at time t-n - t-0,
c                       Ntime=6 => Value at time t-n+alpha - t-0,
c                       Ntime=7 => Value at time t-n - t-0.
c       xi ...........  Nodal transformation matrix at time m
c                       xi(3,3,ElementNode)
c       shp ........... Nodal shape function values.
c       shp1,shp2 ..... Nodal shape function natural derivatives.
c       shx1,shx2 ..... Nodal shape function cartesian derivatives.
c       sg ............ Gauss point coordinates and weights.
c       b ............. Discrete strain-displacement operator.
c       c ............. Resultant stress constitutive matrix.
c       zphG .......... Initial local coordinates
c                       at Gauss points.
c       zph1,zph2 ..... Initial coordinate local derivatives
c                       at Gauss points.
c       zpx1,zpx2 ..... Initial coordinate global derivatives
c                       at Gauss points.
c       zphm .......... Initial coordinate local derivatives
c                       at midside nodes.
c       cphi .......... Current local nodal coordinates.
c       cpx1,cpx2 ..... Current coordinate global derivatives
c                       at Gauss points.
c       cphm .......... Current coordinate local derivatives
c                       at midside nodes.
c       xjinv ......... Local-global coordinate Jacobian matrix.
c       xjs ........... Reference jacobian at Gauss points.
c       xjw ........... Reference jacobian weighted for integration.
c       zdx1,zdx2 ..... Initial director global derivatives
c                       at Gauss points.
c       zdrm .......... Initial local directors
c                       at midside nodes.
c       cdx1,cdx2 ..... Current director global derivatives
c                       at Gauss points.
c       cdrm .......... Current local directors
c                       at midside nodes.
c       ze,ce ......... Initial, current membrane strain measure.
c       zx,cx ......... Initial, current shear strain measure.
c       zr,cr ......... Initial, current bending strain measure.
c       sn ............ Membrane stress.
c       sq ............ Shear stress.
c       sm ............ Bending stress.
c       aa ............ Gauss point accelerations at t-n+a.
c       cmem .......... Membrane constitutive relations coeffiecient
c       cbnd .......... Bending constitutive relations coeffiecient
c       nu ............ Poisson's ratio.
c       rhoa .......... Trans. thickness weighted density multiplyied
c                       by jacobian and weight at Gauss point.
c       rhoi .......... Rot. inertia weighted density multiplyied
c                       by jacobian and weight at Gauss point.
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
      include  'pmod2d.h'
      include  'prstrs.h'
      include  'strnum.h'
      include  'tdata.h'

      include  'comblk.h'

      logical   optmix
      logical   rhs       , lhs        , errck

      integer   i         , j          , k          , l
      integer   i5        , i6         , j5         , j6
      integer   ij        , ik         , lint
      integer   ndf       , ndm        , nst        , isw

      real*8    rhoa      , rhoi       , xjdet      , xlgnrm
      real*8    xj1       , xj2        , xj3
      real*8    xj11      , xj12       , xj21       , xj22
      real*8    cmem      , cbnd
      real*8    dthk      , fac1       , facn

      real*8    f1   (3)  , f2  (3)    , y1  (4)    , y2  (4)
      real*8    d    (*)  , xl  (ndm,*), ul  (ndf,nen,*)
      real*8    s   (nst,*), p   (nst)
      real*8    sg   (3,4), shp (4,4)
      real*8    shp1 (4,4), shp2(4,4)  , shx1(4,4)  , shx2(4,4)
      real*8    fphm (3,4), fphi(3,4)  , ufphm(3,4) , ufphi(3,4)
      real*8    bphm (3,4), bphi(3,4)  , ubphm(3,4) , ubphi(3,4)
      real*8    zphm (3,4), zph1(3,4)  , zph2(3,4)
      real*8    cphm (3,4), cphi(3,4)  , ucphm(3,4) , ucphi(3,4)
      real*8    zpx1 (3,4), zpx2(3,4)  , cpx1(3,4)  , cpx2(3,4)
      real*8    ucpx1(3,4), ucpx2(3,4)
      real*8    ubpx1(3,4), ubpx2(3,4)
      real*8    ufpx1(3,4), ufpx2(3,4)
      real*8    bpx1 (3,4), bpx2(3,4)  , fpx1(3,4)  , fpx2(3,4)
      real*8    zdrm (3,4), zdx1(3,4)  , zdx2(3,4)
      real*8    bdrm (3,4), bdx1(3,4)  , bdx2(3,4)  , ubdrm(3,4)
      real*8    cdrm (3,4), cdx1(3,4)  , cdx2(3,4)  , ucdrm(3,4)
      real*8    fdrm (3,4), fdx1(3,4)  , fdx2(3,4)  , ufdrm(3,4)
      real*8    ubdx1(3,4), ubdx2(3,4)
      real*8    ucdx1(3,4), ucdx2(3,4)
      real*8    ufdx1(3,4), ufdx2(3,4)
      real*8    ze   (3,4), zx  (2,4)  , zr  (3,4)
      real*8    be   (3,4), bx  (2,4)  , br  (3,4)
      real*8    ce   (3,4), cx  (2,4)  , cr  (3,4)
      real*8    fe   (3,4), fx  (2,4)  , fr  (3,4)
      real*8    sn   (3,4), sq  (2,4)  , sm  (3,4)
      real*8    aa   (3)  , xnt  (13)
      real*8    ploc (24) , sloc(24,24), stmp(3,3)
      real*8    c    (8,8), b(8,24,4)  , b1(8,24,4)
      real*8    xjinv(2,2), xjs (4)    , xjw (4)    , xlg (3,3)
      real*8    dir(3,9,7), xi(3,3,9)  ,xj(3,3,9)

      save

c     Set Thickness from Global Input

      if (isw.gt.1) then
        if(d(14).lt.1.d-5) then
          dthk  = 0.25d0 * (thkl(1)+thkl(2)+thkl(3)+thkl(4))
          if(dthk.lt.0.0d0) then
            write(iow,3002) dthk
            call plstop()
          endif
        else
          dthk  = d(14)
        endif
      endif

c     Unused FEAP Task Options (isw=2,7,9,10,11,12,14):

c     Input Material Properties:

      if(isw.eq.1) then

c       Check input data

        errck = .false.
        if(d(1).le.0.0d0) then
          write(iow,3000)
          errck = .true.
        endif
        if(d(14).le.0.0d0) then
          write(iow,3001)
        endif
        if(errck) stop '  Shell material property errors'

c       Initialize Shape Functions.

        l     = 2
        call int2d      ( l       , lint    , sg )
        do i = 1 , lint
          call sh3fshap ( sg (1,i), shp(1,i), shp1(1,i), shp2(1,i)  )
        end do ! i

c       Set default rotational update type

        rotyp = -1

c       Set Number of Printed/Plotted Stresses to 8.

        istv = max(istv,10)

c     Compute Element Dynamic Tangent and Residual Array:

      elseif(isw.ge.3 .and. isw.le.6  .or. isw.eq.8 .or. isw.eq.13) then

c       Set local element director fields

        call sh3finte (xl, ul, xln, ndof, ndm, ndf, nel,  dir, xi, xj)

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
          xj2    =  zph1(3,i)*zph2(1,i) - zph1(1,i)*zph2(3,i)
          xj3    =  zph1(1,i)*zph2(2,i) - zph1(2,i)*zph2(1,i)
          xjs(i) =  sqrt( xj1**2 + xj2**2 + xj3**2 )
          if(isw.eq.4 .or. isw.eq.8) then
            xjw(i) =  1.d0
          else
            xjw(i) =  sg(3,i) * xjs(i)
          endif

c       End Gauss Loop

        end do ! i

c       Do mass computation

        if(isw.eq.5) then
          call sh3fmasg( d, dthk, shp, xjw, s, lint, nst)
          return
        endif

c       Set Flags from Input

        optmix = d(17).gt.1.d0
        lhs    = isw.eq.3
        rhs    = .true.

c       Compute Current Configuration

        fac1 = 1.d0 / theta(3) - 1.d0
        do i = 1 , 4
          do j = 1 , 3
            cphi(j,i)  = xl(j,i)    + ul(j,i,1)
            bphi(j,i)  = cphi(j,i)  - ul(j,i,2)
            fphi(j,i)  = cphi(j,i)  + ul(j,i,2)*fac1
            ucphi(j,i) = ul(j,i,1)
            ubphi(j,i) = ucphi(j,i) - ul(j,i,2)
            ufphi(j,i) = ucphi(j,i) + ul(j,i,2)*fac1
          end do ! j
        end do ! i

c       Compute Midside Interpolations

        call sh3fintr ( xl   , dir(1,1,4) , zphm , zdrm )

        call sh3fintr ( bphi , dir(1,1,1) , bphm , bdrm )
        call sh3fintr ( cphi , dir(1,1,2) , cphm , cdrm )
        call sh3fintr ( fphi , dir(1,1,3) , fphm , fdrm )

        call sh3fintr (ubphi , dir(1,1,5) ,ubphm ,ubdrm )
        call sh3fintr (ucphi , dir(1,1,6) ,ucphm ,ucdrm )
        call sh3fintr (ufphi , dir(1,1,7) ,ufphm ,ufdrm )

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

          call sh3flmda ( xlg(1,3) , xlg  )

          xj11       = xlg(1,1)*zph1(1,i) + xlg(2,1)*zph1(2,i)
     &               + xlg(3,1)*zph1(3,i)
          xj21       = xlg(1,2)*zph1(1,i) + xlg(2,2)*zph1(2,i)
     &               + xlg(3,2)*zph1(3,i)
          xj12       = xlg(1,1)*zph2(1,i) + xlg(2,1)*zph2(2,i)
     &               + xlg(3,1)*zph2(3,i)
          xj22       = xlg(1,2)*zph2(1,i) + xlg(2,2)*zph2(2,i)
     &               + xlg(3,2)*zph2(3,i)

          xjdet      = 1.d0/(xj11*xj22 - xj12*xj21)

          xjinv(1,1) = xj22 * xjdet
          xjinv(1,2) =-xj12 * xjdet
          xjinv(2,1) =-xj21 * xjdet
          xjinv(2,2) = xj11 * xjdet

c         Compute Global Shape Function Derivatives

          do j = 1 , 4
             shx1(j,i) = xjinv(1,1)*shp1(j,i) + xjinv(2,1)*shp2(j,i)
             shx2(j,i) = xjinv(1,2)*shp1(j,i) + xjinv(2,2)*shp2(j,i)
          end do ! j

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

            ucpx1(j,i) = shx1(1,i)*ucphi(j,1)   + shx1(2,i)*ucphi(j,2)
     &                 + shx1(3,i)*ucphi(j,3)   + shx1(4,i)*ucphi(j,4)
            ucpx2(j,i) = shx2(1,i)*ucphi(j,1)   + shx2(2,i)*ucphi(j,2)
     &                 + shx2(3,i)*ucphi(j,3)   + shx2(4,i)*ucphi(j,4)
            ufpx1(j,i) = shx1(1,i)*ufphi(j,1)   + shx1(2,i)*ufphi(j,2)
     &                 + shx1(3,i)*ufphi(j,3)   + shx1(4,i)*ufphi(j,4)
            ufpx2(j,i) = shx2(1,i)*ufphi(j,1)   + shx2(2,i)*ufphi(j,2)
     &                 + shx2(3,i)*ufphi(j,3)   + shx2(4,i)*ufphi(j,4)
            ubpx1(j,i) = shx1(1,i)*ubphi(j,1)   + shx1(2,i)*ubphi(j,2)
     &                 + shx1(3,i)*ubphi(j,3)   + shx1(4,i)*ubphi(j,4)
            ubpx2(j,i) = shx2(1,i)*ubphi(j,1)   + shx2(2,i)*ubphi(j,2)
     &                 + shx2(3,i)*ubphi(j,3)   + shx2(4,i)*ubphi(j,4)

c           Directors

            zdx1(j,i) = shx1(1,i)*dir(j,1,4) + shx1(2,i)*dir(j,2,4)
     &                + shx1(3,i)*dir(j,3,4) + shx1(4,i)*dir(j,4,4)
            zdx2(j,i) = shx2(1,i)*dir(j,1,4) + shx2(2,i)*dir(j,2,4)
     &                + shx2(3,i)*dir(j,3,4) + shx2(4,i)*dir(j,4,4)

            bdx1(j,i) = shx1(1,i)*dir(j,1,1) + shx1(2,i)*dir(j,2,1)
     &                + shx1(3,i)*dir(j,3,1) + shx1(4,i)*dir(j,4,1)
            bdx2(j,i) = shx2(1,i)*dir(j,1,1) + shx2(2,i)*dir(j,2,1)
     &                + shx2(3,i)*dir(j,3,1) + shx2(4,i)*dir(j,4,1)
            cdx1(j,i) = shx1(1,i)*dir(j,1,2) + shx1(2,i)*dir(j,2,2)
     &                + shx1(3,i)*dir(j,3,2) + shx1(4,i)*dir(j,4,2)
            cdx2(j,i) = shx2(1,i)*dir(j,1,2) + shx2(2,i)*dir(j,2,2)
     &                + shx2(3,i)*dir(j,3,2) + shx2(4,i)*dir(j,4,2)
            fdx1(j,i) = shx1(1,i)*dir(j,1,3) + shx1(2,i)*dir(j,2,3)
     &                + shx1(3,i)*dir(j,3,3) + shx1(4,i)*dir(j,4,3)
            fdx2(j,i) = shx2(1,i)*dir(j,1,3) + shx2(2,i)*dir(j,2,3)
     &                + shx2(3,i)*dir(j,3,3) + shx2(4,i)*dir(j,4,3)

            ubdx1(j,i) = shx1(1,i)*dir(j,1,5) + shx1(2,i)*dir(j,2,5)
     &                 + shx1(3,i)*dir(j,3,5) + shx1(4,i)*dir(j,4,5)
            ubdx2(j,i) = shx2(1,i)*dir(j,1,5) + shx2(2,i)*dir(j,2,5)
     &                 + shx2(3,i)*dir(j,3,5) + shx2(4,i)*dir(j,4,5)
            ucdx1(j,i) = shx1(1,i)*dir(j,1,6) + shx1(2,i)*dir(j,2,6)
     &                 + shx1(3,i)*dir(j,3,6) + shx1(4,i)*dir(j,4,6)
            ucdx2(j,i) = shx2(1,i)*dir(j,1,6) + shx2(2,i)*dir(j,2,6)
     &                 + shx2(3,i)*dir(j,3,6) + shx2(4,i)*dir(j,4,6)
            ufdx1(j,i) = shx1(1,i)*dir(j,1,7) + shx1(2,i)*dir(j,2,7)
     &                 + shx1(3,i)*dir(j,3,7) + shx1(4,i)*dir(j,4,7)
            ufdx2(j,i) = shx2(1,i)*dir(j,1,7) + shx2(2,i)*dir(j,2,7)
     &                 + shx2(3,i)*dir(j,3,7) + shx2(4,i)*dir(j,4,7)
          end do ! j

c         Compute Strain Measures

          if (isw.eq.13 .or. .not.fl(9)) then

            call sh3fstrns( shp1(1,i)     , shp2(1,i)   ,
     &                      zphm          , zpx1(1,i)   , zpx2(1,i)   ,
     &                      zdrm          , zdx1(1,i)   , zdx2(1,i)   ,
     &                      ufphm         , ufpx1(1,i)  , ufpx2(1,i)  ,
     &                      ufdrm         , ufdx1(1,i)  , ufdx2(1,i)  ,
     &                      ce  (1,i)     , cx  (1,i)   , cr  (1,i)   )

          else

            call sh3fstrn ( shp1(1,i)     , shp2(1,i)   ,
     &                      zphm          , zpx1(1,i)   , zpx2(1,i)   ,
     &                      zdrm          , zdx1(1,i)   , zdx2(1,i)   ,
     &                      ze  (1,i)     , zx  (1,i)   , zr  (1,i)   )
            call sh3fstrn ( shp1(1,i)     , shp2(1,i)   ,
     &                      fphm          , fpx1(1,i)   , fpx2(1,i)   ,
     &                      fdrm          , fdx1(1,i)   , fdx2(1,i)   ,
     &                      fe  (1,i)     , fx  (1,i)   , fr  (1,i)   )
            call sh3fstrn ( shp1(1,i)     , shp2(1,i)   ,
     &                      bphm          , bpx1(1,i)   , bpx2(1,i)   ,
     &                      bdrm          , bdx1(1,i)   , bdx2(1,i)   ,
     &                      be  (1,i)     , bx  (1,i)   , br  (1,i)   )
            do j = 1 , 3
              ce(j,i) = fac1*fe(j,i) + facn*be(j,i) - ze(j,i)
              cr(j,i) = fac1*fr(j,i) + facn*br(j,i) - zr(j,i)
            end do ! j
            do j = 1 , 2
              cx(j,i) = fac1*fx(j,i) + facn*bx(j,i) - zx(j,i)
            end do ! j
          endif

c         Set Up Strain-Displacement Matrices

          call sh3fbmtx ( shp1(1,i)     , shp2(1,i)   ,
     &                    shx1(1,i)     , shx2(1,i)   ,
     &                    cphm          ,
     &                    cpx1(1,i)     , cpx2(1,i)   ,
     &                    cdrm          ,
     &                    cdx1(1,i)     , cdx2(1,i)   ,
     &                    b(1,1,i)                    )

          if (lhs) then

            call sh3fbmtx ( shp1(1,i)     , shp2(1,i)   ,
     &                      shx1(1,i)     , shx2(1,i)   ,
     &                      fphm          ,
     &                      fpx1(1,i)     , fpx2(1,i)   ,
     &                      fdrm          ,
     &                      fdx1(1,i)     , fdx2(1,i)   ,
     &                      b1(1,1,i)                    )
          endif

c       End Gauss Loop

        end do ! i

c       Zero full (local) stiffness and residuals

        do i = 1,24
          do j = 1,24
            sloc(j,i) = 0.0d0
          end do ! j
          ploc(i) = 0.0d0
        end do ! i

c       General Pian-Sumihara Calculations

        if (optmix) then
          cmem = d(1)*dthk
          cbnd = d(1)*dthk**3/12.d0
          call sh3fpian ( xl   , xjw , sg , cmem , cbnd ,
     &                    d(2) , f1  , f2 , y1 , y2   , xnt  )

c         Pian-Sumihara Membrane Treatment

          call sh3fmbrn ( xjw  , sg   , cmem , d(2) ,
     &                    f1   , f2   , y1  , y2   , xnt  ,
     &                    6    , 24   , ce  , b    , sn   ,
     &                    ploc , sloc )

c         Pian-Sumihara Bending Treatment

          call sh3fbend ( xjw  , sg  , cbnd , d(2) ,
     &                    f1   , f2  , y1 , y2   , xnt  ,
     &                    6    , 24  , cr , b    , sm   ,
     &                    ploc , sloc )
        endif

c       Galerkin and Shear Treatment

        do i = 1,lint

c         Compute Constitutive Relations (weighted by jacobian)

          call sh3fcmtx ( d      , dthk      , xjw(i)    ,
     &                    optmix , zph1(1,i) , zph2(1,i) , c )

c         Compute Stresses (weighted by jacobian)

          call sh3fstrs ( c         ,
     &                    ce(1,i)   , cx(1,i)   , cr(1,i)   ,
     &                    optmix    ,
     &                    sn(1,i)   , sq(1,i)   , sm(1,i)   )

c         Store time history plot data for element

          i6 = 8*(l-1)
          do j = 1,3
            tt(j+i6  ) = sn(j,i)
            tt(j+i6+3) = sm(j,i)
          end do ! j
          tt(7+i6) = sq(1,i)
          tt(8+i6) = sq(2,i)

c         Compute Residual

          call sh3fresd ( b (1,1,i) ,
     &                    sn(1,i)   , sq(1,i)   , sm(1,i)   ,
     &                    optmix    , ploc      )

          if (lhs) then

c           Compute Material Stiffness

            call sh3fstif ( b (1,1,i) , b1(1,1,i) , c         ,
     &                      optmix    , sloc                  )

c           Compute Geometric Stiffness

            if(gflag) then
              call sh3fgeom ( shp1(1,i) , shp2(1,i) ,
     &                        shx1(1,i) , shx2(1,i) ,
     &                        cphm      , cpx1(1,i) , cpx2(1,i) ,
     &                        dir(1,1,2),
     &                        sn  (1,i) , sq  (1,i) , sm  (1,i) ,
     &                        sloc                              )
            endif
          endif

c       End Gauss Loop

        end do ! i

c       Compute Transient Terms

        if (fl(9) .and. noi.ne.6 .and. noi.gt.0) then

          do i = 1,lint

c           Interpolate Accelerations

            do j = 1 , 3

              aa(j) = shp(1,i)*ul(j,1,5) + shp(2,i)*ul(j,2,5)
     &              + shp(3,i)*ul(j,3,5) + shp(4,i)*ul(j,4,5)

            end do ! j

c           Compute Weighted Density/Inertia

            rhoa = d(4) * d(14) * xjw(i)
            rhoi = d(8) * rhoa * d(14)**2/12.d0

c           Compute Transient Residual and Tangent

            call sh3ftran ( dt,  shp(1,i)         , aa               ,
     &                      rots(1,1,1)           , rots(1,1,2)      ,
     &                      rvel(1,1,1)           , rvel(1,1,2)      ,
     &                      racc(1,1,1)           , racc(1,1,2)      ,
     &                      dir(1,1,2)            , dir(1,1,4)       ,
     &                      rhoa         , rhoi   , ploc      , sloc )

c         End Gauss Loop

          end do ! i

        endif

c       Stress and energy processes

        if (isw.eq.4 .or. isw.eq.8)  then
          call sh3fstre ( xl,shp, ce,cx,cr,sn,sq,sm,xjw, lint,ndm,isw )
          rhs = .false.
        elseif (isw.eq.13) then
          call sh3fener ( d,dthk, ce,cr,cx,sn,sq,sm,
     &                    ul,shp,fphi,xjw, ndf,lint )
          rhs = .false.
        elseif(lhs) then

c         Scale tangent stiffness for total increment du_t-n+1

          do j = 1,24
            do i = 1,24
              sloc(i,j) = fac1*sloc(i,j)
            end do ! i
          end do ! j

c         Transform (Reduce) Stiffness/Residual - First Node Loop

          i5 = 0
          i6 = 0

c         First Node Loop

          do i = 1 , 4

            j5 = 0
            j6 = 0

c           Second Node Loop

            do j = 1 , 4

c             Stiffness Translational Terms

              do l = 1 , 3
                do k = 1 , 3
                  s(i5+k,j5+l) = sloc(i6+k,j6+l)
                end do ! k
              end do ! l

c             Stiffness Off-diagonal Terms

              do l = 1 , ndof(j) - 3
                do k = 1 , 3
                  s(i5+k,j5+l+3) = sloc(i6+k,j6+4)*xj(1,l,j)
     &                           + sloc(i6+k,j6+5)*xj(2,l,j)
     &                           + sloc(i6+k,j6+6)*xj(3,l,j)
                end do ! k
              end do ! l
              do l = 1 , ndof(i) - 3
                do k = 1 , 3
                  s(i5+l+3,j5+k) = sloc(i6+4,j6+k)*xi(1,l,i)
     &                           + sloc(i6+5,j6+k)*xi(2,l,i)
     &                           + sloc(i6+6,j6+k)*xi(3,l,i)
                end do ! k
              end do ! l

c             Stiffness Rotational Terms

              do l = 1 , ndof(j) - 3
                do k = 1 , 3
                  stmp(k,l) = sloc(i6+k+3,j6+4)*xj(1,l,j)
     &                      + sloc(i6+k+3,j6+5)*xj(2,l,j)
     &                      + sloc(i6+k+3,j6+6)*xj(3,l,j)
                end do ! k
              end do ! l
              do l = 1 , ndof(j) - 3
                do k = 1 , ndof(i) - 3
                  s(i5+k+3,j5+l+3) = stmp(1,l)*xi(1,k,i)
     &                             + stmp(2,l)*xi(2,k,i)
     &                             + stmp(3,l)*xi(3,k,i)
                end do ! k
              end do ! l

              j5 = j5 + ndf
              j6 = j6 + 6

c           End Second Node Loop

            end do ! j

            i5 = i5 + ndf
            i6 = i6 + 6

c         End First Node Loop

          end do ! i

        endif

c       Node Loop

        if(rhs) then

          i5 = 0
          i6 = 0

          do i = 1 , 4

c           Residual Terms

            do j = 1 , 3
              p(i5+j) = ploc(i6+j)
              ij      = i6 + 3 + j
              do k = 1 , ndof(i) - 3
                ik = i5 + 3 + k
                p(ik) = p(ik) + ploc(ij)*xi(j,k,i)
              end do ! k
            end do ! j

            i5 = i5 + ndf
            i6 = i6 + 6

c         End Node Loop

          end do ! i

c         Add geometric stiffness term for 6 dof nodes

          if(lhs .and. ndf.eq.6 .and. gflag)then
            call sh3fgeo6 (lint , shp1 , shp2 , shx1 , shx2 ,
     &                    cphm , cpx1 , cpx2 ,  dir (1,1,2) ,
     &                    sq  , sm    , s                   )
          end if

c         Modify residual and stiffnes for pressure & body loading

          call sh3fsurf(d,dthk,xl,ul,xjw, ndf,ndm,nst, p,s)

        endif
      endif

c     Formats:

3000  format(' *ERROR* Modulus is zero or negative for shell')

3001  format(' *WARNING* Thickness zero: Compute from nodal input.')

3002  format(' *ERROR* Thickness is zero or negative for shell')

      end
