c$Id:$
      subroutine framf3c(d,ul,xl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Comment all rotational inertia terms             03/09/2008
c       2. Remove call to pltln -- not required for plots   14/10/2008
c       3. Use nn2,nn3 as increment to nh2                  17/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Static finite deformation, geometrically exact rod

c     Implementated by: Ignacio Romero (iromero@ce.berkeley.edu)
c                       November 2002

c     Based on an idea first described in Crisfield and Jelenic:
c     Crisfield, M. A. and Jeleni\'c}[1998], "Objectivity of strain
c     measures in geometrically exact 3D beam theory and its finite
c     element implementation", Proceedings of the Royal Society of
c     London, Series A, 455, 1983, 1125--1147

c     Only SVK constitutive model. Other material models can be added
c     in modelsvq.
c     Only for static problems.

c     Usage notes:
c         use rotational updates type 8 (exp map), included in file.
c         Spatial dimension      : 3
c         Number of dof per node : 6
c         Number of nodes/element: ONLY 2

c     Notes:
c      - mass matrix needs to be corrected
c      - The current element is both objective and path independent (for
c        elastic problems)

c     Functions in this file:
c        framf3c           : main file

c        bmatcj            : linearized strain matrix computation
c        checkcj           : check element geometry for errors
c        commonscj         : common blocks definitions
c        enercj            : energy and momenta calculation
c        forcecj           : body forces, follower and conservative
c        hatcj             : skew symmetric matrix associated to vector
c        hmatcj            : matrix appearing in linearization of expmap
c        initializecj      : compute initial rotations
c        interpmatrix      : interpolation matrix of increment variation
c        kgeocj            : geometric tangent
c        kmatcj            : material tangent
c        masscj            : mass matrix computation
c        modelcj           : stress and material elasticities
c        momentcj          : spatial section momentum
c        positioncj        : centroid position computation
c        pushfcj           : push forward of stress and tangent
c        slerp             : spherical interpolation
c        statcj            : static residual and tangent
c        straincj          : calculation of strain measures
c        strenergycj       : stored energy function
c        tangcj            : tangent, residual and stress computation
c        tmatcj            : compute matrix 't', in linearization
c        updrotationcj     : update of rotation tensors
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      implicit   none

      include   'eldata.h' !n, nel
      include   'erotas.h' !rotyp
      include   'hdata.h'  !nh
      include   'iofile.h' !iow

      logical    symflg
      integer    ndf, ndm, nst, isw
      real*8     d(*),ul(ndf,*),xl(ndm,*),s(nst,nst),p(nst)

      integer    mint, gint

      save

      if(isw.eq.0) then
        if(ior.lt.0) then
          write(*,*) '   3d nonlinear frame, Crisfield & Jelenic'
        endif

c     Input element parameters, done in inmate.f

      else if(isw.eq.1) then

c       Set amount of memory in history variables per element

        if (nel .eq. 2) then
          gint = 1
          mint = 2
        else
          gint = 2
          mint = 3
        end if
        if(nint(d(182)).gt.0) then      ! Nodal quadrature
          gint = mint
        endif
        nh1 = nh1*gint

c       In these history variables we store initial rotation tensor

        nh3 = 4*(gint+mint)

c       This element (will) work with two rotational updates:
c          Cayley transform(?)
c          Exponential map (-2).
c       We pick first one only for conserving solution, else, second
c       this has to be changed, so element recognizes when solution is
c       conserving

        rotyp =-2

c       flag to warn about symmetric tangent

        symflg = .true.

c     Check for zero length and for bad reference nodes

      else if (isw .eq. 2) then

        call checkcj(n, ndm, nel, xl)

      else

c       Initialization of common blocks

        call commonscj(d)

c       Computation of tangent stiffness and residual vector

        if (isw.eq.3 .or. isw.eq.6 .or. isw.eq.4 .or. isw.eq.8) then

          call tangcj(d, ul, xl, s, p, ndf, ndm, nst, isw, symflg)

c       Compute element mass matrices

        elseif(isw.eq.5) then

          call masscj(d, xl , nel , ndf , ndm , nst , s , p)

c       Update of database just before start of a new time step
c       Nothing to be done for C-J rod

        elseif(isw.eq.12) then

c       Computes and prints angular momentum and kinetic, potential
c       and total energies

        elseif(isw.eq.13) then

c         energy computations disabled
c         call enercj(nel, ndf, ndm, d , xl, ul)

c       Initialize database values. Puts rotations in undeformed
c       configuration inside tn and tn+1 position of history variables

        elseif(isw.eq.14) then

          call initializecj(d, xl, ndm, nel, n)

        endif

      endif ! isw

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: bmatcj

c     Purpose:    computes linearized strain matrix B

c     Input:      nel   : number of nodes in element
c                 shp   : shape functions and derivatives
c                 phipr : derivative w/r to S of centroid curve
c                 qmat  : interpolation matrix in linearized rotation
c                 dmat  : matrix in linearization of curvature

c     Output:     bmat  : matrix B
c                 brmat : matrix Br, for right multiplication in Kmat
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine bmatcj(nel, shp, phipr, qmat, dmat, bmat, brmat)

      implicit none

      integer  nel
      real*8   shp(2,*), phipr(3), qmat(3,3,*), dmat(3,3)
      real*8   bmat(6,6,*), brmat(6,6,*)

      integer  a, i,j
      real*8   hp(3,3)

c     Compute standard matrix Bmat

      do a = 1,nel

        do i = 1,6
          do j = 1,6
            bmat(j,i,a) = 0.0d0
          end do ! j
          bmat(i,i,a) = shp(1,a)
        end do ! i

        bmat(1,5,a) = -shp(2,a)*phipr(3)
        bmat(1,6,a) =  shp(2,a)*phipr(2)
        bmat(2,4,a) =  shp(2,a)*phipr(3)
        bmat(2,6,a) = -shp(2,a)*phipr(1)
        bmat(3,4,a) = -shp(2,a)*phipr(2)
        bmat(3,5,a) =  shp(2,a)*phipr(1)

      end do ! a

c     Compute modified matrix Brmat

      do a = 1,nel

        do i = 1,6
          do j = 1,6
            brmat(i,i,a) = 0.0d0
          end do ! j
        end do ! i

c       du-du block

        do i = 1,3
          brmat(i,i,a) = shp(1,a)
        end do ! i

c       du-dtheta block

        call hatcj( phipr, hp)
        do j = 1,3
          do i = 1,3
            brmat(i,j+3,a) = hp(i,1)*qmat(1,j,a)
     &                     + hp(i,2)*qmat(2,j,a)
     &                     + hp(i,3)*qmat(3,j,a)
          end do ! i
        end do ! j

c       dthetha-dtheta block

        do j = 1,3
          do i = 1,3
            brmat(i+3,j+3,a) = dmat(i,j)*shp(1,a)
          end do ! i
        end do ! j

      end do ! a

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: checkcj

c     Purpose:    checks element for input errors in geometry

c     Input:      n   : node number
c                 ndm : spatial dimension
c                 nel : nodes in element
c                 xl  : reference nodal coordinates

c     Output:     only output messages

c     Notes :     only for 2 noded elements
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine checkcj(n, ndm, nel, xl)

      implicit   none

      include   'iofile.h'
      include   'refnd.h'

      integer    n, ndm, nel
      real*8     xl(ndm, *)

      integer    k
      real*8     inc(3) , norm  , length
      real*8     t1(3)  , t2(3) , t3(3)

c     First check: zero element length

      if (nel .eq. 2) then

        do k = 1,3
          inc(k) = xl(k,2) - xl(k,1)
        end do ! k
        norm = inc(1)*inc(1) + inc(2)*inc(2) + inc(3)*inc(3)

        if ( norm .lt. 1d-5) then
          write (iow, 2000) n
          if (ior .lt. 0) then
            write(*, 2000) n
          endif
        end if

c       Second check, test for bad reference node/vector

        length = 1.d0/sqrt(norm)
        do k = 1,3
          t3(k) = inc(k)*length
        end do ! k

        if    (lref.eq.1) then ! REFE NODE
          do k = 1,3
            t2(k) = refx(k) - xl(k,1)
          end do ! k
        elseif(lref.eq.2) then ! REFE VECTor
          do k = 1,3
            t2(k) = refx(k)
          end do ! k
        elseif(lref.eq.3) then ! REFE POLAr
          t2(1)   = 0.5d0*(xl(1,1)+ xl(1,2))
          t2(2)   = 0.5d0*(xl(2,1)+ xl(2,2))
          t2(3)   = atan2(t2(2),t2(1))
          length  = sqrt(t2(1)*t2(1) + t2(2)*t2(2))
          t2(1)   = length*cos(t2(3))
          t2(2)   = length*sin(t2(3))
          t2(3)   = 0.0d0
        endif

        call vecp(t2,t3,t1)
        norm = t1(1)*t1(1) + t1(2)*t1(2) + t1(3)*t1(3)

c       Check for bad REFErence data.
        if ( sqrt(norm) .lt. 1.d-8) then

          write(iow,2100) n
          if (ior .lt. 0)  write(*,2100) n

        end if
      end if

 2000 format(//' Warning: Element ', i4, ' has zero length')
 2100 format(//' Warning: bad REFErence node/vector for element', i4)

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: commonscj

c     Purpose: Sets variables in /bmcj1/,/bmcj2/ from element
c              properties

c     Input:   d: material data

c     Output:  everything goes in common blocks

c     Notes:   At each integration points total rotation stored. That
c              is, rotation tensor that transform inertial frame to
c              section frame.  Observe that there is also posibility
c              to store incremental rotation tensor that only accounts
c              for transformation between section frame at t=0 and
c              current one.
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine commonscj(d)

      implicit   none

      include   'bm3fd.h'
      include   'cdata.h'
      include   'eldata.h'
      include   'eqsym.h'
      include   'refnd.h'
      include   'sdata.h'

      integer    i
      real*8     d(*)

c     Set assembly vector

      do i = 1,6
        is(i   ) = i
        is(i+ 6) = i + ndf
        is(i+12) = i + ndf + ndf
      end do ! i

c     Select integration rule for stiffness (gint) & mass (mint) term.
c     Use reduced integration for stiffness to avoid locking.

      if(nel.eq.2) then
        gint   = 1
        mint   = 2
      else
        gint   = 2
        mint   = 3
      endif

c     Select Gauss-Legendre or Gauss-Lobatto integration type.
c     To select Lobatto use 'QUAD,,n' with n<0 in material defintion.

      if(nint(d(6)).lt.0) then
        lobatt = 1
      else
        lobatt = 0
      endif

c     Set flags for symmetrize tangent or not

      if(neqs.eq.neq) then
        isym = 1
      else
        isym = 0
      endif

c     Computation of mct0, number of history terms for element.
c     per Gauss point in stiffnes: 4 (rotation quaternion)_n + 3
c                                     spatial curvature omega_n
c     per Gauss point in mass    : 4 (rotation quaternion)_n
c      mct0 = 7*gint+4*mint
c     no history variables
      mct0  = 0

c     Material inertia matrix and weighted area

      miner(1,1) = d(4)*d(33)*d(8)
      miner(1,2) = d(4)*d(35)*d(8)
      miner(1,3) = 0.0d0
      miner(2,1) = miner(2,1)
      miner(2,2) = d(4)*d(34)*d(8)
      miner(2,3) = 0.0d0
      miner(3,1) = 0.0d0
      miner(3,2) = 0.0d0
      miner(3,3) = miner(1,1) + miner(2,2)
      arho       = d(4)*d(32)

c     Set distributed loading type.
c     Use 'GFOL' or 'LFOL' in material definition
c     'LOCA' -> d(67) = 2;  conservative load defined in rod local axis
c     'GFOL' -> d(67) = 3;  follower load in global axis
c     'LFOL' -> d(67) = 4;  follower load in local axis
c     else   -> d(67) = 1;  conservative load in global axis

      iforc  = int(d(67))

c     External loading, multiplied by proportional value

      eforce(1) = d(11)*dm
      eforce(2) = d(12)*dm
      eforce(3) = d(13)*dm

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: enercj

c     Purpose:    computation of energy and momentum

c     Input:      nel: nodes in element
c                 ndf: degrees of freedom per node
c                 d  : feap material data
c                 xl : nodal reference positions
c                 ul : degrees of freedom vector

c     Output:     In file P*.ene:

c                 epl(1-3) : linear momentum
c                 epl(4-6) : angular momentum
c                 epl(7)   : Kinetic energy
c                 epl(8)   : Potential energy
c                 epl(9)   : External energy
c                 epl(10)  : Total energy             (in pmacr2.f)
c                 epl(11)  : Angular momentum norm    (in pmacr2.f)

c     Notes:      for isw=13, independently of integrator,
c                 all variables (u,v,a) in ul() are at time tn+1
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine enercj(nel , ndf , ndm , d , xl , ul)

      implicit   none

      include   'bm3fd.h'
      include   'cdata.h'  !nen
      include   'comblk.h' !hr
      include   'erotas.h' !rvel,
      include   'ddata.h'  !noi
      include   'hdata.h'  !nh
      include   'ptdat6.h' !epl

      integer    nel         , ndf       , ndm
      real*8     xl(ndm,*)   , ul(ndf,*) , d(*)

      integer    nn2         , ip        , kk       , a
      integer    nns         , nhs       , nh       , isw
      real*8     sg(2,5)     , shp(2,3)  , xjac
      real*8     dx          , w
      real*8     xnode1(3,3) , omega1(3)
      real*8     phi1(3)     , phipr1(3)
      real*8     qrot1(4)
      real*8     rotmat1(3,3), strain(6)

c     Nodal positions at time tn+1

      do a = 1,nel
        do kk = 1,3
          xnode1(kk,a) = xl(kk,a) + ul(kk,a)
        end do ! kk
      end do ! a

c     Computation of Gauss points and weights for stiffness integration

      if(nint(d(182)).gt.0) then
        call int1dn(gint, sg)
      else
        call int1d(gint, sg)
      endif

c     Loop on each Gauss point

      nn2 = 0
      nns = mct0
      nh  = nint(d(15))  ! Number of history variables/quad. point
      nhs = nint(d(149)) ! Number of history variables/section
      do ip = 1, gint

c       Recover data from history variables
c       (rotation and curvature at tn+1)

        do kk = 1,4
          qrot1(kk)  = hr(nh2+nn2-1+kk)
        end do ! kk
        do kk = 1,3
          omega1(kk) = hr(nh2+nn2+3+kk)
        end do ! kk

c       Shape functions and derivatives. Arclength element

        call shp1d(sg(1,ip), xl, shp, ndm, nel, xjac)
        dx = sg(2,ip)*xjac

c       Computation of position vectors & derivative with respect to S

        call positioncj(xnode1, nel, shp, phi1, phipr1)

c       Computation of rotation matrix at time tn+1

        call quamat(qrot1 , rotmat1)

c       Compute convected strains.

c       Compute stored energy

        call strenergycj(d , hr(nh1+nns), hr(nh2+nns), nh,
     &                   strain, w, isw)

c       accumulate stored enegy, integrated

        epl(8) = epl(8) + w*dx

        nn2 = nn2 + 7
        nns = nns + nhs

      end do ! ip: gauss loop

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: forcecj

c     Purpose:    compute distributed load

c     Input:      rotmatr: initial rotation tensor
c                 rotmat: rotation tensor
c                 eforce: body force, as input by user
c                 iforc : type of load.
c                         (1) : conservative load given in global axis
c                         (2) : conservative loald in rod local axis
c                         (3) : follower load in global axis
c                         (4) : follower load in rod local axis

c     Output:     bodyf : rotated body force
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine forcecj(rotmatr, rotmat, eforce, iforc, bodyf)

      implicit   none

      integer    iforc, k
      real*8     rotmatr(3,3), rotmat(3,3), eforce(3), bodyf(3), tmp(3)

c     conservative load in global axis

      if (iforc .eq. 1) then

        do k = 1,3
          bodyf(k) = eforce(k)
        end do ! k

c     conservative load in local axis = Lambda_o*eforce

      else if (iforc .eq. 2) then

        do k = 1,3
          bodyf(k) = rotmatr(k,1)*eforce(1)
     &             + rotmatr(k,2)*eforce(2)
     &             + rotmatr(k,3)*eforce(3)
        end do ! k

c     follower loads in global axis = Lambda*Lambda_o^t*eforce

      else if (iforc .eq. 3) then

        do k = 1,3
          tmp(k) = rotmatr(1,k)*eforce(1)
     &           + rotmatr(2,k)*eforce(2)
     &           + rotmatr(3,k)*eforce(3)
        end do ! k

        do k = 1,3
          bodyf(k) = rotmat(k,1)*tmp(1)
     &             + rotmat(k,2)*tmp(2)
     &             + rotmat(k,3)*tmp(3)
        end do ! k

c     follower loads in local axis

      else if (iforc .eq. 4) then

        do k = 1,3
          bodyf(k) = rotmat(k,1)*eforce(1)
     &             + rotmat(k,2)*eforce(2)
     &             + rotmat(k,3)*eforce(3)
        end do ! k

      end if

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: hatcj

c     Purpose:    compute skew symmetric matrix associated to R3 vector

c     Input:      v  :(3)    vector

c     Output:     m  :(3,3)  matrix
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine hatcj(v, m)

      implicit   none

      real*8     v(3), m(3,3)

      m(1,1) =  0.0d0
      m(1,2) = -v(3)
      m(1,3) =  v(2)
      m(2,1) =  v(3)
      m(2,2) =  0.0d0
      m(2,3) = -v(1)
      m(3,1) = -v(2)
      m(3,2) =  v(1)
      m(3,3) =  0.0d0

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: hmatcj

c     Purpose:    Compute matrix h appearing when linearizing exp map

c     Input:      psi : rotation vector

c     Output:     h

c     Notes:      h = nn + sin(b)/b (eye-nn) + (1-cos(b))/n hatn,
c                 b = norm(psi),
c                 n = psi/b,
c                 nn= n otimes n
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine hmatcj(psi, h)

      implicit   none

      real*8     psi(3), h(3,3)

      integer    i , j
      real*8     b, s, c, z
      real*8     hn(3,3)

      b = psi(1)*psi(1) + psi(2)*psi(2) + psi(3)*psi(3)

      if ( b .lt. 1.0d-12) then
        s = (1.0d0 - 0.05d0*b)/6.0d0
        c =  0.5d0 - b*(1.d0 - b/30.d0)/24.d0
        b =  1.0d0 - b*s
      else
        b = sqrt(b)
        z = sin(b)/b
        c = (1.d0 - cos(b))/b**2
        s = (1-z)/b**2
        b = z
      end if

      call hatcj(psi, hn)
      do i = 1,3
        do j = 1,3
          h(j,i) = s*psi(j)*psi(i) + c*hn(j,i)
        end do ! j
        h(i,i) = h(i,i) + b
      end do ! i

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: initializecj

c     Purpose:    Initializes database at start of whole process

c     Input:      d   : feap material data
c                 xl  : reference nodal coordinates
c                 ndm : spatial dimension
c                 nel : nodes in element
c                 n   : node number

c     Output:     goes in common block

c     Notes:      History variables initialized in following order:
c                   -- first all Gauss points for stiffness calculation
c                   -- then, all Gauss points for mass calculation
c                 For each G.p. for stiffness,  following data is stored
c                   -- a quaternion, with rotation
c                   -- spatial curvature, a (3) vector
c                 For each G.p. for mass:
c                   -- a quaternion, with rotatin

c                 Important: initialization routine commonscj must be
c                            called before this one
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine initializecj(d, xl, ndm, nel, n)

      implicit   none

      include   'bm3fd.h'
      include   'comblk.h'
      include   'erotas.h'
      include   'hdata.h'
      include   'refnd.h'
      include   'iofile.h' !iow

      integer    ndm , nel , n
      real*8     d(*), xl(ndm, *)

      integer    ip, k, nn3
      integer    j1,k1,l1
      real*8     t1(3),t2(3),t3(3),length,sg(2,5)
      real*8     xlbda0(3,3),shp(2,3),xjac,aa
      real*8     sgl(2,5)

      save

c     External node, and coordinates
c     In material definition, one must give one of following
c     REFE, NODE , x1, x2, x3 ... and coordinates of node  t1, or
c     REFE, VECT,  v1, v2, v3 ... and components of vector t1
c     REFE, POLA,  v1, v2, v3 ... and components of vector t1

      lref    = nint(d(96))
      refx(1) = d(97)
      refx(2) = d(98)
      refx(3) = d(99)

c     Checking of straight element

      if(nel.eq.3) then

        if((xl(1,3)-xl(1,1)).eq.0.d0 .and.
     &     (xl(1,3)-xl(1,2)).ne.0.d0) go to 101

        if((xl(2,3)-xl(2,1)).eq.0.d0 .and.
     &     (xl(2,3)-xl(2,2)).ne.0.d0) go to 101

        if((xl(3,3)-xl(3,1)).eq.0.d0 .and.
     &     (xl(3,3)-xl(3,2)).ne.0.d0) go to 101

        do j1 = 1,3
          if((xl(j1,3)-xl(j1,1)).ne.0.d0) then
            t1(j1) = (xl(j1,3)-xl(j1,2))/(xl(j1,3)-xl(j1,1))
          else
            if(j1.ne.1) t1(j1) = t1(1)
          end if
        end do ! j1
        if((xl(1,3)-xl(1,1)).ne.0.d0) then
          aa    = t1(1)
        else
          aa    = t1(2)
          t1(1) = aa
        end if
        do j1 = 1,3
          t1(j1) = (t1(j1)-aa)/aa
        end do ! j1
        if(t1(2).gt.1.d-03.or.t1(3).gt.1.d-03) go to 101
      end if

c     Straight element. Computation of unique rotation matrix.
c     This is only part for two noded elements

      do j1 = 1,3
        t3(j1) = xl(j1,2) - xl(j1,1)
      end do ! j1
      length = 1.d0/sqrt(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3))
      do j1 = 1,3
        t3(j1) = t3(j1)*length
      end do ! j1

      if    (lref.eq.1) then
        do j1 = 1,3
          t2(j1) = refx(j1) - xl(j1,1)
        end do ! j1
      elseif(lref.eq.2) then
        do j1 = 1,3
          t2(j1) = refx(j1)
        end do ! j1
      elseif(lref.eq.3) then ! REFE POLAr
        t2(1)   = 0.5d0*(xl(1,1)+ xl(1,2))
        t2(2)   = 0.5d0*(xl(2,1)+ xl(2,2))
        t2(3)   = atan2(t2(2),t2(1))
        length  = sqrt(t2(1)*t2(1) + t2(2)*t2(2))
        t2(1)   = length*cos(t2(3))
        t2(2)   = length*sin(t2(3))
        t2(3)   = 0.0d0
      endif

      call vecp(t2,t3,t1)
      length = sqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))

c     Check for bad REFErence data.

      if ( length .lt. 1.d-8) then

         write(iow,*)
         write(iow,*) 'Warning: bad REFErence node/vector for element',n
         write(iow,*) 'Invalid computations. Correct mesh.'

         if (ior .lt. 0) then
           write(*,*) 'Warning: bad REFErence node/vector for element',n
         endif

      else

         length = 1.d0/length
         do j1 = 1,3
           t1(j1) = t1(j1)*length
         end do ! j1
         call vecp(t3,t1,t2)

      end if

c     Form rotation matrix and extract quaternion

      do k1 = 1,3
        xlbda0(k1,1) = t1(k1)
        xlbda0(k1,2) = t2(k1)
        xlbda0(k1,3) = t3(k1)
      end do ! k1

c     pass rot tensor to vector and store one per int. point in rotqu0
c     in quaternion format

      call quatex(xlbda0,rotqu0(1,1))

      do j1 = 2,gint+mint
        do k1 = 1,4
          rotqu0(k1,j1) = rotqu0(k1,1)
        end do ! k1
      end do ! j1

      goto 200

c     Curved element

c     Gauss points for stifness

  101 if(nint(d(182)).gt.0) then
        call int1dn(gint, sg)
      else
        call int1d(gint, sg)
      endif

      do j1 = 1,gint
        call shp1d(sg(1,j1), xl, shp, ndm, nel, xjac)
        do k1 = 1,3
          t3(k1) = 0.0d0
        end do ! k1

c       The third vector has direction of tangent to line of centroids

        do l1 = 1,nel
          do k1 = 1,3
            t3(k1) = t3(k1) + shp(1,l1)*xl(k1,l1)
          end do ! k1
        end do ! l1
        length = 1.d0/sqrt(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3))
        do k1 = 1,3
          t3(k1) = t3(k1)*length
        end do ! k1

        if    (lref.eq.1) then
          do k1 = 1,3
            t2(k1) = refx(k1) - xl(k1,1)
          end do ! k1
        elseif(lref.eq.2) then
          do k1 = 1,3
            t2(k1) = refx(k1)
          end do ! k1
        elseif(lref.eq.3) then ! REFE POLAr
          t2(1)   = 0.5d0*(xl(1,1)+ xl(1,2))
          t2(2)   = 0.5d0*(xl(2,1)+ xl(2,2))
          t2(3)   = atan2(t2(2),t2(1))
          length  = sqrt(t2(1)*t2(1) + t2(2)*t2(2))
          t2(1)   = length*cos(t2(3))
          t2(2)   = length*sin(t2(3))
          t2(3)   = 0.0d0
        endif
        call vecp(t2,t3,t1)
        length = 1.d0/sqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
        do k1 = 1,3
          t1(k1) = t1(k1)*length
        end do ! k1
        call vecp(t3,t1,t2)

c       Transfer to rotqu0

        do k1 = 1,3
          xlbda0(k1,1) = t1(k1)
          xlbda0(k1,2) = t2(k1)
          xlbda0(k1,3) = t3(k1)
        end do ! k1
        call quatex(xlbda0,rotqu0(1,j1))
      end do ! j1

c     Gauss points for mass

      if(lobatt.eq.0) then
        call int1d (mint,sgl)
      else
        call int1dn(mint,sgl)
      end if

      do j1 = 1,mint
        call shp1d(sgl(1,j1),xl,shp,ndm,nel,xjac)
        do k1 = 1,3
          t3(k1) = 0.0d0
        end do ! k1

c       The third vector has direction of tangent to line of centroids

        do l1 = 1,nel
          do k1 = 1,3
            t3(k1) = t3(k1) + shp(1,l1)*xl(k1,l1)
          end do ! k1
        end do ! l1
        length = 1.d0/sqrt(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3))
        do k1 = 1,3
          t3(k1) = t3(k1)*length
        end do ! k1

        if    (lref.eq.1) then
          do k1 = 1,3
            t2(k1) = refx(k1) - xl(k1,1)
          end do ! k1
        elseif(lref.eq.2) then
          do k1 = 1,3
            t2(k1) = refx(k1)
          end do ! k1
        elseif(lref.eq.3) then ! REFE POLAr
          t2(1)   = 0.5d0*(xl(1,1)+ xl(1,2))
          t2(2)   = 0.5d0*(xl(2,1)+ xl(2,2))
          t2(3)   = atan2(t2(2),t2(1))
          length  = sqrt(t2(1)*t2(1) + t2(2)*t2(2))
          t2(1)   = length*cos(t2(3))
          t2(2)   = length*sin(t2(3))
          t2(3)   = 0.0d0
        endif
        call vecp(t2,t3,t1)
        length = 1.d0/sqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
        do k1 = 1,3
          t1(k1) = t1(k1)*length
        end do ! k1
        call vecp(t3,t1,t2)

c       Transfer to rotqu0

        do k1 = 1,3
          xlbda0(k1,1) = t1(k1)
          xlbda0(k1,2) = t2(k1)
          xlbda0(k1,3) = t3(k1)
        end do ! k1
        call quatex(xlbda0,rotqu0(1,j1+gint))
      end do ! j1

c     This rotation matrices are stored in nh3, they never change
c     because we will need to know initial rotation

 200  nn3 = 0
      do ip = 1,gint+mint
         do k = 1,4
            hr(nh3+nn3-1+k) = rotqu0(k,ip)
         end do ! ik
         nn3 = nn3 + 4
      end do ! ip

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: interpmatrix

c     Purpose:    Compute interpolation matrices needed for consistent
c                 incremental rotation

c     Input:      qnode :  rotation (quaternion) at nodes
c                 qi :  interpolated rotation at gauss point
c                 qref: reference rotation
c                 sl :  interpolation factor in [0,1]

c     Output:     qmat: 2 interpolation matrices for rotation increments
c                 dmat: matrix for curvature increment

c     Notes:      qbar    = Lambda_1*h(t*psi)*t(psi)*Lambda1^t,
c                 qmat(1) = eye-q2
c                 qmat(2) = t*qbar,
c                 dmat    = Lambda*Lambdaref^t*t(psi)*Lambda1^t
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine interpmatrix(qnode, qi, qref, sl, qmat, dmat)

      implicit   none

      real*8     qnode(4,2), qi(4), qref(4), sl
      real*8     qmat(3,3,2) , dmat(3,3)

      real*8     psi(3), psit(3) , lambda1(3,3), h(3,3), t(3,3)
      real*8     qrel(4), qtmp(4)
      real*8     tmp1(3,3), tmp2(3,3)
      integer    a, i, j

c     Compute rotation angle from lambda1 to lambda2

      call qutrmul(qnode(1,1), qnode(1,2), qrel, 1)
      call quarot(qrel , psi)

c     Compute rotation angle from lambda1 to rotation at gauss point

      do a = 1,3
        psit(a) = psi(a) * sl
      end do ! a

c     Compute factors

      call quamat(qnode(1,1), lambda1)
      call hmatcj(psit, h)
      call tmatcj(psi , t)

c     Multiply factors to get interpolating matrices

      do j = 1,3
        do i = 1,3
          tmp1(i,j) = h(i,1)*t(1,j) + h(i,2)*t(2,j) + h(i,3)*t(3,j)
        end do ! i
      end do ! j

      do j = 1,3
        do i = 1,3
          tmp2(i,j) = lambda1(i,1)*tmp1(1,j)
     &              + lambda1(i,2)*tmp1(2,j)
     &              + lambda1(i,3)*tmp1(3,j)
        end do ! i
      end do ! j
      do j = 1,3
        do i = 1,3
          qmat(i,j,2) = sl*(tmp2(i,1)*lambda1(j,1)
     &                    + tmp2(i,2)*lambda1(j,2)
     &                    + tmp2(i,3)*lambda1(j,3))
          qmat(i,j,1) = -qmat(i,j,2)
        end do ! k
        qmat(j,j,1) = qmat(j,j,1) + 1.d0
      end do ! j

c     Computation of second matrix D

      call qutrmul(qi  , qref , qtmp, 2)
      call quamat(qtmp, tmp1)

      do j = 1,3
        do i = 1,3
          tmp2(i,j) = tmp1(i,1)*t(1,j)
     &              + tmp1(i,2)*t(2,j)
     &              + tmp1(i,3)*t(3,j)
        end do ! i
      end do ! j

      do j = 1,3
        do i = 1,3
          dmat(i,j) = tmp2(i,1)*lambda1(j,1)
     &              + tmp2(i,2)*lambda1(j,2)
     &              + tmp2(i,3)*lambda1(j,3)
        end do ! i
      end do ! j

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: kgeocj

c     Purpose:    compute geometric part of tangent stiffness

c     Input:      cons: true if conserving scheme
c                 nel : nodes on element
c                 shp : shape functions and derivatives
c                 phipr : derivative w/r to S of centroid curve
c                 qmat: interpolation matrices for rotation increments
c                 sstress: spatial stresses [n m]

c     Output:     kgeo : tangent
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine kgeocj( nel, shp, phipr, qmat, sstress, kgeo)

      implicit   none

      integer    nel
      real*8     shp(2,*), phipr(3), qmat(3,3,2) , sstress(6)
      real*8     kgeo(18,18)

      integer    a, b, aa, bb , i, j
      real*8     hatn(3,3) , hatm(3,3)
      real*8     f(3,3) , x

      call hatcj( sstress(1), hatn )
      call hatcj( sstress(4), hatm )

      x = sstress(1)*phipr(1)+sstress(2)*phipr(2)+sstress(3)*phipr(3)
      do i = 1,3
         do j = 1,3
           f(j,i) = sstress(j)*phipr(i)
         end do ! j
         f(i,i) = f(i,i) - x
      end do ! i

      aa = 0
      do a = 1,nel
        bb = 0
        do b = 1,nel

c         Block 1-1

          do j = 1,3
            do i = 1,3
              kgeo(aa+i,bb+j) = 0.0d0
            end do ! i
          end do ! j

c         Block 1-2

          do j = 1,3
            do i = 1,3
              kgeo(aa+i,bb+3+j) =  -shp(1,a)*(hatn(i,1)*qmat(1,j,b)
     &                                      + hatn(i,2)*qmat(2,j,b)
     &                                      + hatn(i,3)*qmat(3,j,b))
            end do ! i
          end do ! j

c         Block 2-1

          do j = 1,3
            do i = 1,3
              kgeo(aa+3+i,bb+j) =   shp(2,a)*shp(1,b)*hatn(i,j)
            end do ! i
          end do ! j

c         Block 2-2

          do j = 1,3
            do i = 1,3
              kgeo(aa+3+i,bb+3+j) = - shp(1,a)*(hatm(i,1)*qmat(1,j,b)
     &                                        + hatm(i,2)*qmat(2,j,b)
     &                                        + hatm(i,3)*qmat(3,j,b))
     &                              + shp(2,a)*(   f(i,1)*qmat(1,j,b)
     &                                        +    f(i,2)*qmat(2,j,b)
     &                                        +    f(i,3)*qmat(3,j,b))
            end do ! i
          end do ! j

          bb = bb+6

        end do ! b
        aa = aa+6
      end do ! a

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: kmatcj

c     Purpose:    compute material part of tangent stiffness

c     Input:
c                 nel : nodes on element
c                 bmat:  linearized strain matrix at equilibrium point
c                 brmat: b matrix for right multiplication
c                 scc : spatial material elasticity

c     Output:     kmat : tangent

c     Notes:      For conserving scheme, Bmat is (Bmat0+Bmat1)/2
c                 For rest is Bmat(Phi_n+alpha)
c                 K = bmat^t cc brmat
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine kmatcj(nel , bmat , brmat , scc, kmat)

      implicit   none

      integer    nel
      real*8     bmat(6,*), brmat(6,*), scc(6,6), kmat(18,18)

      integer    j,k
      real*8     btc(6,18)

c     Multiply Bt C and store its transpose

      do j = 1,6*nel
        do k = 1,6
          btc(k,j) = bmat(1,j)*scc(1,k)
     &             + bmat(2,j)*scc(2,k)
     &             + bmat(3,j)*scc(3,k)
     &             + bmat(4,j)*scc(4,k)
     &             + bmat(5,j)*scc(5,k)
     &             + bmat(6,j)*scc(6,k)
        end do ! k
      end do ! j

c     Complete multiplication. Not asumed symmetric

      do k = 1,6*nel
        do j = 1,6*nel
          kmat(j,k) = btc(1,j)*brmat(1,k)
     &              + btc(2,j)*brmat(2,k)
     &              + btc(3,j)*brmat(3,k)
     &              + btc(4,j)*brmat(4,k)
     &              + btc(5,j)*brmat(5,k)
     &              + btc(6,j)*brmat(6,k)
        end do ! j
      end do ! k

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: masscj

c     Purpose:    computation of mass matrix

c     Input:      d  : material parameters
c                 xl : nodal reference positions
c                 nel: nodes in element
c                 ndm: spatial dimension
c                 nst: dimension of mass matrix (ndm*nel)

c     Output:     s: mass matrix
c                 p: lumped mass matrix

c     Notes:      mass matrix is not constant for this model
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine masscj(d, xl , nel , ndf , ndm , nst , s , p)

      implicit   none

      include   'bm3fd.h'
      include   'comblk.h'
      include   'cdata.h'
      include   'ddata.h'  !noi, theta
      include   'eltran.h'
      include   'erotas.h'
      include   'hdata.h'
      include   'tdata.h'

      integer    ndf ,  ndm       , nst        , nel
      real*8     d(*),  xl(ndm,*) , s(nst,nst) , p(nst)

      logical    cons
      integer    a, b, i, k1,l1, ip
      integer    nn1

      real*8     dx         , xjac
      real*8     sg(2,5)    , shp(2,3)
      real*8     qrot0(4)   , qrot(4)     , qrot1(4)
      real*8     thts(3)    , tht1(3)
      real*8     rotmat1(3,3)
      real*8     wn1(3)     , pin1(3)     , in1(3,3)
      real*8     hpin1(3,3)
      real*8     tmatrix(3,3)
      real*8     tmp        , tmp1(3,3)   , tmp2(3,3)

      cons = noi .eq. 5

c     Computation of Gauss points and weights for stiffness integration

      if(nint(d(182)).gt.0) then
        call int1dn(mint, sg)
      else
        call int1d(mint, sg)
      endif

c     Loop on each Gauss point

      nn1   = 7*gint

      do ip = 1,mint

c       Recover data from history variables (rotation_n)

        do i = 1,4
          qrot0(i)  = hr(nh1 + nn1-1+i)
        end do ! i

c       Shape functions and derivatives. Arclength element

        call shp1d(sg(1,ip), xl, shp, ndm, nel, xjac)
        dx = sg(2,ip)*xjac

c       Computation of spatial rotation vector from Lambda_n to
c       Lambda_n+alpha

        do i = 1,3
          thts(i) = 0.0d0
          do a = 1,nel
            thts(i) = thts(i) + theta(3)*rots(i,a,2)*shp(2,a)
            tht1(i) = tht1(i) +          rots(i,a,2)*shp(2,a)
          end do ! a
        end do ! i

c       Computation of rotation_n+1, rotation_n+a at Gauss point
c       (quaternion for

        call updrotationcj(cons, qrot0 , tht1 , qrot1)
        call updrotationcj(cons, qrot0 , thts , qrot )

c       Computation of rotation tensor (matrices) at time tn and tn+1

        call quamat(qrot1 , rotmat1)

c       Angular velocities

        do i = 1,3
          wn1(i) = 0.d0
          do a = 1,nel
            wn1(i) = wn1(i) + shp(2,a) * rvel(i,a,2)
          end do ! a
        end do ! i

c       Calculations for rotational block

        call hatcj(pin1 , hpin1)
        call tmatcj(thts, tmatrix)

        do a = 1,3
          do b = 1,3
            tmp1(a,b) = tmatrix(b,a)
          end do ! b
        end do ! a

c       multiply spatial inertia (symmetric) by tmp1

        do b = 1,3
          do a = 1,3
            tmp2(a,b) = in1(1,a)*tmp1(1,b)
     &                + in1(2,a)*tmp1(2,b)
     &                + in1(3,a)*tmp1(3,b)
          end do ! a
        end do ! b

        l1 = 0
        do a = 1,nel

          k1 = 0
          do b = 1,nel

c           translation contribution

            tmp          = shp(2,a)*shp(2,b)*dx*ctan(3)*arho

            s(l1+1,k1+1) = s(l1+1,k1+1) + tmp
            s(l1+2,k1+2) = s(l1+2,k1+2) + tmp
            s(l1+3,k1+3) = s(l1+3,k1+3) + tmp

c           Rotation contribution

c           tmp          = shp(2,a)*shp(2,b)*dx*ctan(3)

c           do kk = 1,3
c             do i = 1,3
c               s(l1+i+3,k1+kk+3) = s(l1+i+3,k1+kk+3) + 1.0d0
c     &                           + tmp2(i,kk)*ctan(3)
c     &                           + hpin1(kk,i)/(dt*ctan(3))
c             end do ! i
c           end do ! kk

            k1 = k1 + ndf
          end do ! b
          l1 = l1 + ndf
        end do ! a

c       Lumped mass matrix, just add up columns

        do a = 1,nst
          do i = 1,nst
            p(a) = p(a) + s(i,a)
          end do ! i
        end do ! a

c       advance history variable pointers

        nn1 = nn1 + 4

      end do ! ip: gauss loop

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: modelcj

c     Purpose:    Computes convected stresses and elasticity matrix

c     Input:      d       : feap's material data
c                 strain  : (6) convected strain [Gamma , Omega ]
c                 hn(*)   : History variables at t_n
c                 nh      : Number of history variables/point
c                 isw     : Switch value from element

c     Output:     stress  : (6) convected stress [N , M]
c                 cc      : (6,6) elasticity tensor
c                 h1(*)   : History variables at t_n+1

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine modelcj(d, hn, h1, nh, strain, stress, cc, isw)

      implicit   none

      integer    i, j, nh, isw

      real*8     d(*),hn(*),h1(*)
      real*8     strain(6), stress(6), cc(6,6)

      real*8     e, nu, a, ix, iy, jz, xy, kx, ky, g

      save

c     General stress/strain relation for general section beams.

      if(d(101).gt.0.0d0) then

        call bm3res(d,hn,h1,nh,strain, stress,cc, isw)

c     Elasticity tensor. Constant for SVK material model

      else
        e  = d(1)
        nu = d(2)
        a  = d(32)
        ix = d(33)
        iy = d(34)
        xy = d(35)
        jz = d(36)
        kx = d(37)
        ky = d(38)
        g  = 0.5d0*e/(1.d0+nu)

        do j = 1,6
          do i = 1,6
            cc(i,j) = 0.0d0
          end do ! i
        end do ! j
        cc(1,1) = g*a*kx
        cc(2,2) = g*a*ky
        cc(3,3) = e*a
        cc(4,4) = e*ix
        cc(5,5) = e*iy
        cc(6,6) = g*jz

        cc(4,5) = e*xy
        cc(5,4) = e*xy

c       Stress resultants. Formula valid only for SVK model

        do i = 1,6
          stress(i) = cc(i,1)*strain(1) + cc(i,2)*strain(2)
     &              + cc(i,3)*strain(3) + cc(i,4)*strain(4)
     &              + cc(i,5)*strain(5) + cc(i,6)*strain(6)
        end do ! i

      endif

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: positioncj

c     Purpose:    Computes current position of one point on line of
c                 centroid (Phi) & derivative with respect to S (Phi,s)

c     Input:      xnode : (3,nel) nodal position vectors
c                 nel   : number of nodes in element
c                 shp   : (2,nel) derivative & shape function

c     Output:     phi   : (3) position vector at Gauss point
c                 phipr : (3) derivative of phi w/r reference arclength

c     Notes:      -
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine positioncj(xnode , nel , shp , phi , phipr)

      implicit   none

      integer    nel, a, k
      real*8     xnode(3,nel), shp(2,nel), phi(3), phipr(3)

      do k = 1,3
        phi(k)   = 0.0d0
        phipr(k) = 0.0d0
        do a = 1,nel
          phi(k)   = phi(k)   + shp(2,a)*xnode(k,a)
          phipr(k) = phipr(k) + shp(1,a)*xnode(k,a)
        end do ! a
      end do ! k

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: pushfcj

c     Purpose:    Compute algorithmic spatial stresses and material
c                 elasticity by pushing forward convective tensors

c     Input:
c                 rotmat : (3,3) rotation matrix at equilibrium
c                 stress : (6)   convected stresses [N , M]
c                 cc     : (6,6) convected elasticities

c     Output:     sstress: (6)   spatial stresses [n , m]
c                 scc    : (6,6) spatial elasticities
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine pushfcj(rotmat, stress, cc, sstress , scc)

      implicit   none

      include   'ddata.h' !theta

      real*8     rotmat(3,3), stress(6) , cc(6,6)
      real*8     sstress(6) , scc(6,6)

      integer    j         , k
      real*8     tmp(6,6)

c     Push forward stresses

      do k = 1,3
         sstress(k  ) = rotmat(k,1)*stress(1)
     &                + rotmat(k,2)*stress(2)
     &                + rotmat(k,3)*stress(3)
         sstress(k+3) = rotmat(k,1)*stress(4)
     &                + rotmat(k,2)*stress(5)
     &                + rotmat(k,3)*stress(6)
      end do ! k

c     Push forward elasticity tensor.
c     Compute transformation matrix for left matrix

c     Multiply rotmat *cc

      do k = 1,6
        do j = 1,3
          tmp(j  ,k) = rotmat(j,1)*cc(1,k)
     &               + rotmat(j,2)*cc(2,k)
     &               + rotmat(j,3)*cc(3,k)
          tmp(j+3,k) = rotmat(j,1)*cc(4,k)
     &               + rotmat(j,2)*cc(5,k)
     &               + rotmat(j,3)*cc(6,k)
        end do ! j
      end do ! k

c     Multiply tmp * rotmat^t

      do k = 1,3
        do j = 1,6
          scc(j,k  ) = tmp(j,1)*rotmat(k,1)
     &               + tmp(j,2)*rotmat(k,2)
     &               + tmp(j,3)*rotmat(k,3)
          scc(j,k+3) = tmp(j,4)*rotmat(k,1)
     &               + tmp(j,5)*rotmat(k,2)
     &               + tmp(j,6)*rotmat(k,3)
        end do ! j
      end do ! k

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: slerp

c     Purpose:    spherical interpolation of rotations

c     Input:      q1: rotation in quaternion format (vector, scalar)
c                 q2:
c                 t : interpolation factor in [0,1]

c     Output:     qi: interpolated rotation in quaternion format

c     Notes:      q1 = q1 [ q1^* q2 ]^t
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine slerp(q1, q2, t, qi)

      implicit   none
      real*8     q1(4), q2(4), t
      real*8     qi(4)

      real*8     qrel(4) , qpower(4)
      real*8     angle, un
      integer    a

c     Compute relative rotation from q1 to q2

      call qutrmul(q1, q2, qrel, 1)

c     Raise qrel to power t

      angle = acos(qrel(4))
      un    = sqrt(qrel(1)*qrel(1) + qrel(2)*qrel(2) + qrel(3)*qrel(3))
      if ( un .gt. 1e-8) then
        qpower(4) = cos(angle*t)
        do a = 1,3
          qpower(a) = qrel(a) * sin(angle*t) / un
        end do ! a
      else
        do a = 1,3
          qpower(a) = 0.0d0
        end do ! a
        qpower(4) = 1.d0
      end if

c     Last operation for slerp

      call quamul(q1 , qpower , qi)
      call quanrm(qi)

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: statcj

c     Purpose:    Compute static part of residual and tangent

c     Input:      d     : feap's material data
c                 ul    : feap's dof vector
c                 xl    : feap's nodal reference positions
c                 ndf   : # of degrees of freedom
c                 ndm   : spatial dimension
c                 nst   : dimension of s and p
c                 isw   : task
c                 xnode0: nodal positions at tn
c                 xnode : nodal positions at tn+alpha (equil point)
c                 xnode1: nodal positions at tn+1

c     Output:     s     : tangent matrix
c                 p     : residual
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine statcj(d, ul, xl, s, p, ndf ,ndm, nst, isw, xnode)

      implicit   none

      include   'bm3fd.h'
      include   'cdata.h'
      include   'comblk.h'
      include   'counts.h'
      include   'ddata.h'  !noi
      include   'debugs.h'
      include   'eldata.h'
      include   'elplot.h' !tt
      include   'eltran.h'
      include   'erotas.h' ! xln, rots
      include   'fdata.h'  ! fl
      include   'hdata.h'  ! nh3
      include   'iofile.h' ! iow, ior
      include   'tdata.h'

      integer    ndf      , ndm       , nst       , isw
      real*8     d(*)     , xl(ndm,*) , ul(ndf,*) , s(nst,nst), p(nst)
      real*8     xnode(3,*)

      integer    a, b, k1,l1, kk, ip, nn3, nh, nhs, nns

      real*8     dx        , xjac        , length
      real*8     sg(2,5)   , shp(2,3)
      real*8     phi(3)    , phipr(3)    , pp(6,2)
      real*8     sl
      real*8     qnode(4,2)
      real*8     qgauss(4) , qref(4)     , qslerp(4)
      real*8     qmat(3,3,2), dmat(3,3)
      real*8     rotmatr(3,3), rotmat(3,3)
      real*8     bodyf(3)
      real*8     strain(6)
      real*8     stress(6) , sstress(6)
      real*8     cc(6,6)   , scc(6,6)
      real*8     bmat(6,6,3), brmat(6,6,3)
      real*8     kmat(18,18), kgeo(18,18)
      real*8     tmp

      save

c     Element reference length

      length = sqrt((xl(1,2)-xl(1,1))*(xl(1,2)-xl(1,1))
     &            + (xl(2,2)-xl(2,1))*(xl(2,2)-xl(2,1))
     &            + (xl(3,2)-xl(3,1))*(xl(3,2)-xl(3,1)))

c     Nodal rotation quaternions

      do a = 1,4
        qnode(a,1) = xln(a,1,3)
        qnode(a,2) = xln(a,2,3)
      end do ! a

c     Computation of Gauss points and weights for stiffness integration

      if(nint(d(182)).gt.0) then
        call int1dn(gint, sg)
      else
        call int1d(gint, sg)
      endif

c     Loop on each Gauss point

      nn3 = 0
      nns = mct0
      nh  = nint(d(15))  ! Number of history variables/quad. point
      nhs = nint(d(149)) ! Number of history variables/section
      do ip = 1, gint

c       Recover data from history variables (rotation_ref)

        do a = 1,4
          qref(a)  = hr(nh3+nn3-1+a)
        end do ! a

c       Shape functions and derivatives. Arclength element

        call shp1d(sg(1,ip), xl, shp, ndm, nel, xjac)
        dx = sg(2,ip)*xjac

c       Computation of position vectors and derivative with respect to S

        call positioncj(xnode , nel, shp, phi , phipr )

c       Value of arclength ratio S/L at gauss point = N_2(xi)

        sl = shp(2,2)

c       Compute rotation tensor at Gauss point: spherical interpolation

        call slerp  (qnode(1,1) , qnode(1,2), sl , qslerp)
        call quamul (qslerp     , qref  , qgauss)
        call quamat (qgauss     , rotmat)

c       Interpolation matrices for incremental rotation

        call interpmatrix(qnode, qgauss, qref, sl, qmat, dmat)

c       Compute convected strains and linearized strain matrix

        call straincj(phipr, rotmat , qnode, qref, length, strain)
        call bmatcj(nel, shp, phipr, qmat, dmat, bmat, brmat)

c       Compute convected stresses & material tangent: push to spatial

        call modelcj(d, hr(nh1+nns), hr(nh2+nns), nh,
     &               strain, stress, cc, isw)
        call pushfcj(rotmat , stress, cc, sstress, scc)

c       Compute contribution from distributed forces (conservative or
c       follower) and contribution to tangent

        call quamat(qref, rotmatr)
        call forcecj(rotmatr, rotmat, eforce, iforc, bodyf)

c       Static contribution to residual

        k1 = 0
        do a = 1,nel

c         Stress divergence terms

          do kk = 1,6
            p(k1+kk) = p(k1+kk) - (bmat(1,kk,a)*sstress(1)
     &                           + bmat(2,kk,a)*sstress(2)
     &                           + bmat(3,kk,a)*sstress(3)
     &                           + bmat(4,kk,a)*sstress(4)
     &                           + bmat(5,kk,a)*sstress(5)
     &                           + bmat(6,kk,a)*sstress(6))*dx
          end do ! kk

c         Body forces (no body couples) contribution

          do kk = 1,3
            p(k1+kk) = p(k1+kk) + shp(2,a)*bodyf(kk)*dx
          end do ! kk

          k1  = k1 + ndf
        end do ! a

c       Computation of tangent

        if (isw .eq. 3) then

          dx = dx * ctan(1)

c         Material part of tangent

          call kmatcj(nel, bmat, brmat, scc, kmat)

c         Geometric part

          call kgeocj(nel, shp, phipr, qmat, sstress, kgeo)

c         Assemble both contributions in element tangent

          do k1 = 1,nel*6
            do l1 = 1,nel*6
              s(is(l1),is(k1)) = s(is(l1),is(k1))
     &                         + (kgeo(l1,k1) + kmat(l1,k1))*dx
            end do ! l1
          end do ! k1

c         Add contributions from follower loads if needed. Only
c         a skew block in du dtheta part

          if (iforc .ge. 3) then

            l1 = 0
            do a = 1,nel

              k1 = 3
              do b = 1,nel

                tmp          = shp(2,a)*shp(2,b)*dx

                s(l1+1,k1+2) = s(l1+1,k1+2) - tmp*bodyf(3)
                s(l1+2,k1+1) = s(l1+2,k1+1) + tmp*bodyf(3)

                s(l1+1,k1+3) = s(l1+1,k1+3) + tmp*bodyf(2)
                s(l1+3,k1+1) = s(l1+3,k1+1) - tmp*bodyf(2)

                s(l1+2,k1+3) = s(l1+2,k1+3) - tmp*bodyf(1)
                s(l1+3,k1+2) = s(l1+3,k1+2) + tmp*bodyf(1)

                k1 = k1 + ndf
              end do ! b
              l1 = l1 + ndf
            end do ! a

          end if

        end if

c       Output of initial coordinates, and configuration

        if (isw .eq. 4) then
          write(iow,2001) n,ip,phi,rotmat
          write(iow,2002) (strain(k1),k1=1,6),(stress(k1),k1=1,6)

          if(ior.lt.0) then
            write(*,2001) n,ip,phi,rotmat
            write(*,2002)(strain(k1),k1=1,6),(stress(k1),k1=1,6)
          endif

        end if

        nn3 = nn3 + 4
        nns = nns + nhs

      end do ! ip: Gauss loop

c     Compute rotational body force residual and tangent

      if(d(4).gt.0.0d0 .and. d(65).gt.0.0d0) then
        call fbody3w(d(4)*d(32),d(65),xl,ul, p,s, isw.eq.3)
      endif

c     If we use only 1 integration point for reduced integration we
c     can save strains and stresses at middle of element.
c     Store strains and stresses for 'tplot,stress,nn,1:12' output, with
c     nn number of element

      do a = 1,6
        tt(a)   = strain(a)
        tt(a+6) = stress(a)
      end do ! a

c     Project nodal end forces to nodes

      if(isw.eq.8) then
        do kk = 1,6
          pp(kk,1) =  p(kk)
          pp(kk,2) = -p(kk+ndf)
        end do ! kk
        call frcn3d(pp,p,s)
      endif

c     Formats

 2001 format(5x,'Element No.',i7,': Point No.',i2/
     &       5x,'Current Coordinates'/5x,1p,3e15.5/
     &       5x,'Current Intrinsic Frame'/5x,3(1p,3e15.5/5x))

 2002 format(5x,'Material Strains'/
     &       5x,'    Gamma : ',1p,3e15.5/
     &       5x,'    Omega : ',1p,3e15.5/
     &       5x,'Material Stresses'/
     &       5x,'    N     : ',1p,3e15.5/
     &       5x,'    M     : ',1p,3e15.5)
      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: straincj

c     Purpose:    Computes convected strain tensors

c     Input:      phipr   : derivative w/r to S of centroid curve
c                 rotmat  : (3,3) rotation matrix
c                 qnode   : rotations at nodes (quaternion)
c                 qref    : undeformed configuration rotn (quaterion)
c                 length  : length of element in ref. configuration

c     Output:     strain  : (6) convected strain [Gamma Omega]
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine straincj(phipr, rotmat, qnode, qref, length,  strain)

      implicit   none

      integer    i
      real*8     phipr(3), rotmat(3,3), qnode(4,2), qref(4), length
      real*8     strain(6)

      real*8     qrel(4), psi(3)

c     Convected strain Gamma
      do i = 1,3
        strain(i) = rotmat(1,i)*phipr(1)
     &            + rotmat(2,i)*phipr(2)
     &            + rotmat(3,i)*phipr(3)
      end do ! i
      strain(3) = strain(3) - 1.d0

c     Convected strain Omega

      call qutrmul(qnode(1,1), qnode(1,2), qrel, 1)
      call quarot (qrel, psi)
      call qutrvec(qref, psi, strain(4))
      do i = 4,6
        strain(i) = strain(i)/length
      end do ! i

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: strenergycj

c     Purpose:    Computes stored energy function

c     Input:      d       : feap's material data
c                 strain  : convected strains
c                     1-3 : strain Gamma
c                     4-6 : strain Omega

c     Output:     w       : stored energy function

c     Notes:      Only valid for quadratic potentials, uses modelcj
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine strenergycj(d , hn, h1, nh, strain , w, isw)
      implicit   none

      integer    nh   , isw
      real*8     d(*) , hn(*), h1(*), strain(6)
      real*8     w

      real*8     stress(6) , cc(6,6)

      call modelcj(d, hn, h1,nh, strain, stress, cc, isw)
      w = 0.5d0*(strain(1)*stress(1) + strain(2)*stress(2)
     &         + strain(3)*stress(3) + strain(4)*stress(4)
     &         + strain(5)*stress(5) + strain(6)*stress(6))

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: tangcj

c     Purpose:    Computes tangent stiffness and residual forces

c     Input:      d     : feap's material data
c                 ul    : feap's dof vector
c                 xl    : feap's nodal reference positions
c                 ndf   : # of degrees of freedom
c                 ndm   : spatial dimension
c                 nst   : dimension of s and p

c     Output:     s     : tangent matrix
c                 p     : residual
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine tangcj(d, ul, xl, s, p, ndf, ndm, nst, isw, symflg)

      implicit   none

      include   'bm3fd.h'
      include   'cdata.h'
      include   'ddata.h'
      include   'eldata.h'
      include   'eltran.h'
      include   'erotas.h' !xln
      include   'fdata.h'  !fl
      include   'tdata.h'

      logical    symflg
      integer    ndf         , ndm        , nst , isw
      integer    a, j, k

      real*8     d(*)        , ul(ndf,*)  , xl(ndm,nel)
      real*8     s(nst,nst)  , p(nst)
      real*8     xnode0(3,3) , xnode(3,3) , xnode1(3,3)
      real*8     fac

      save

c     Nodal positions at t=n+alpha, tn and tn+1

      fac = 1.d0/theta(3) - 1.d0
      do a = 1,nel
        do k = 1,3
          xnode(k,a) = xl(k,a) + ul(k,a)

          if(fl(9)) then
            xnode0(k,a) = xnode(k,a) - ul(k,a+nen)
            xnode1(k,a) = xnode(k,a) + ul(k,a+nen)*fac
          endif

        end do ! k
      end do ! a

c     Compute static  tangent stiffness and residual vector

      call statcj(d, ul, xl, s, p, ndf, ndm , nst, isw, xnode)

c     Symmetrization of stiffness matrix

      if(isym.eq.1) then

        if (symflg) then
          write(*,*) '  Symmetrizing tangent'
          symflg = .false.
        endif

        do j = 1,nst-1
          do k = j+1,nst
            s(j,k) = (s(j,k) + s(k,j))*0.5d0
            s(k,j) =  s(j,k)
          end do ! k
        end do ! j
      end if

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: tmatcj

c     Purpose:    compute matrix t in linearization of rotation vector

c     Input:      thts : rotation vector

c     Output:     tmatrix

c     Notes:      this matrix is different if updates are done with
c                 exponential map or Cayley transform.
c                 The one here is for exp map
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine tmatcj(thts, tmatrix)

      implicit   none

      include   'erotas.h' !rotyp

      real*8     thts(3), tmatrix(3,3)

      integer    i   , j
      real*8     norm, f, c, s, z, htheta(3,3)

      norm = thts(1)*thts(1) + thts(2)*thts(2) + thts(3)*thts(3)

      if ( norm .lt. 1.d-12) then

        z =  0.25d0*norm
        s = (1.d0 + 0.2d0*z*(2.d0 + 17.d0*z/21.d0))/3.d0
        f =  1.d0/(1.d0 + s*z)
        c = -s*f

      else

        norm = sqrt(norm)
        f    = 0.5d0*norm/(tan(0.5d0*norm))
        c    = (1.d0 - f)/(norm*norm)

      endif

      call hatcj(thts, htheta)
      do i = 1,3
        do j=1,3
          tmatrix(j,i) =c*thts(j)*thts(i) - 0.5d0*htheta(j,i)
        end do ! j
        tmatrix(i,i) = tmatrix(i,i) + f
      end do ! i

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: updrotationcj

c     Purpose:    Update rotation matrix given a certain spatial
c                 incremental rotation vector,
c                 using exponential map -> rot2 = exp[tht] rot1
c                 or Cayley transform   -> rot2 = cay[tht] rot1

c     Input:      cons: true if conserving method, uses Cayley
c                 rot   : (4) rotation tensor in quaternion format
c                 tht   : (3) rotation vector

c     Output:     rot2  : (4) updated rotation tensor, quaternion format
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine updrotationcj(cons , qrot0 , tht , qrot1)

      implicit   none

      logical    cons
      real*8     qrot0(4), tht(3), qrot1(4), qtht(4)

      save

c     Computation of quaternions associated to tht

      if (cons) then
        call rqcay( tht, qtht)
      else
        call rotqua( tht, qtht)
      end if

c     Computation of quaternions of new rotation matrix

      call quamul(qtht, qrot0, qrot1)
      call quanrm(qrot1)

      end
