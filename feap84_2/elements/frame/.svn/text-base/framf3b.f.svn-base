c$Id:$
      subroutine framf3b(d,ul,xl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove 'ul' array from subprogram masssvq        13/12/2006
c       2. Comment all rotational inertia terms             03/09/2008
c       3. Remove call to pltln -- not required for plots   14/10/2008
c       4. Use 'nn1,nn3' as increments to nh1,nh2,nh3       17/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Static finite deformation, geometrically exact rod

c     Implementated by: Ignacio Romero (iromero@ce.berkeley.edu)
c                       November 2002

c     Based on Simo & Vu Quoc element of:
c         Simo, J. C. and Vu-Quoc, L. [1986], "A three-dimensional
c         finite-strain rod model. Part II: Computational aspects",
c         Computer Methods in Applied Mechanics and Engineering,
c         58, 79--116

c     Only for static problems.

c     Usage notes:
c         use rotational updates type 8 (exp map), included in file.
c         Spatial dimension      : 3
c         Number of dof per node : 6
c         Number of nodes/element: 2 or 3

c     Notes:
c      - 3 noded elements may not fully work yet. Not tested
c      - 3 noded elements don't draw well
c      - mass matrix needs to be corrected
c      - The paper by Simo & Vuquoc has one serious inconsistency.
c        Throughtout paper, rotation vector from time tn to tn+1
c        is interpolated (denoted \theta).
c        If this is so, linearization is wrong.
c        It can be checked that Lin[ Lambda_n+1 ] \neq delta theta
c        Lambda_n+1 and, among other things, material tangent not
c        symmetric. If theta means rather rotation vector from k-th
c        iteration to k+1-th iteration of rotation, then linearization
c        is ok. This last approach is one followed in current work.
c      - As noted in Crisfield and Jelenic [1998], current element is
c        neither objective nor path invariant

c     Functions in this file:
c     framf3b            : main file

c     bmatsvq            : linearized strain matrix computation
c     checksvq           : check element geometry for errors
c     commonssvq         : common blocks definitions
c     enersvq            : energy and momenta calculation
c     expmapsvq          : exponential map of a vector
c     forcesvq           : body forces, follower and conservative
c     hatsvq             : skew symmetric matrix associate to vector
c     initializesvq      : compute initial rotations
c     kgeosvq            : geometric tangent
c     kmatsvq            : material tangent
c     masssvq            : mass matrix computation
c     modelsvq           : stress and material elasticities
c     positionsvq        : centroid position computation
c     pushfsvq           : push forward of stress and tangent
c     statsvq            : static residual and tangent
c     strainsvq          : calculation of strain measures
c     strenergysvq       : stored energy function
c     tangsvq            : tangent, residual and stress computation
c     tmatsvq            : computes matrix 'T' appears in linearization
c     updomegasvq        : update of spatial curvature
c     updrotationsvq     : update of rotation tensors
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      implicit   none

      include   'eldata.h'       !n, nel
      include   'erotas.h'       !rotyp
      include   'hdata.h'        !nh
      include   'iofile.h'       !iow

      logical    symflg
      integer    ndf, ndm, nst, isw
      real*8     d(*),ul(ndf,*),xl(ndm,*),s(nst,nst),p(nst)

      integer    gint, mint

      save

      if(isw.eq.0) then
        if(ior.lt.0) then
          write(*,*) '   FRM3FD: 3-d Frame Finite Defm (2/3 node)'
        endif

c     Input element parameters, done in inmate.f

      else if(isw.eq.1) then

c       Set amount of memory in history variables per element
c       in stiffness gauss points 7 doubles (4 for rotation,
c       3 for stress) in mass gauss points 4 doubles (for rotation)

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

        nh1 = (nh1 + 7)*gint + 4*mint  ! Change nh1 for history section
                                       ! to history for whole element.

c       In these history variables we store initial rotation tensor

        nh3 = nh3 + 4*(gint+mint)

c       element works with 'prot02'

        rotyp = -2

c       flag to warn about symmetric tangent

        symflg = .true.

c     Check for zero length and for bad reference nodes

      else if (isw .eq. 2) then

        call checksvq(n, ndm, nel, xl)

c     Other options on isw

      else

c       Initialization of common blocks

        call commonssvq(d)

c       Computation of tangent stiffness and residual vector

        if (isw.eq.3 .or. isw.eq.6 .or. isw.eq.4 .or. isw.eq.8) then

          call tangsvq(d, ul, xl, s, p, ndf, ndm, nst, isw, symflg)

c       Compute element mass matrices

        elseif(isw.eq.5) then

          call masssvq(xl , nel , ndf , ndm , nst , s , p)

c       Update of database just before start of a new time step
c       Nothing to be done for S.VQ rod

        elseif(isw.eq.12) then

c       Computes and prints angular momentum and kinetic, potential
c       and total energies

        elseif(isw.eq.13) then

          call enersvq(nel, ndf, ndm, d , xl, ul, isw)

c       Initialize database values. Puts rotations in undeformed
c       configuration inside tn and tn+1 position of history variables

        elseif(isw.eq.14) then

          call initializesvq(d, xl, ndm, nel, n)

        endif

      endif ! isw tests

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: adjugatesvq

c     Purpose:    computes adjugate of a 3x3 matrix

c     Input:      m     : 3x3 matrix

c     Output:     madj  : adjugate of m,
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine adjugatesvq(m, madj)

      implicit   none

      real*8     m(3,3), madj(3,3)

      madj(1,1) = m(2,2)*m(3,3) - m(2,3)*m(3,2)
      madj(1,2) = m(2,3)*m(3,1) - m(2,1)*m(3,3)
      madj(1,3) = m(2,1)*m(3,2) - m(2,2)*m(3,1)

      madj(2,1) = m(3,2)*m(1,3) - m(1,2)*m(3,3)
      madj(2,2) = m(1,1)*m(3,3) - m(1,3)*m(3,1)
      madj(2,3) = m(1,2)*m(3,1) - m(1,1)*m(3,2)

      madj(3,1) = m(1,2)*m(2,3) - m(2,2)*m(1,3)
      madj(3,2) = m(2,1)*m(1,3) - m(1,1)*m(2,3)
      madj(3,3) = m(2,2)*m(1,1) - m(1,2)*m(2,1)

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: averagesvq

c     Purpose:    computes convex average of two vectors

c     Input:      size  : length of vectors
c     alpha : see below
c     v1    : first vector
c     v2    : second vector

c     Output:     w     : average = (1-alpha)*v1 + alpha*v2
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine averagesvq(size, alpha, v1, v2, w)

      implicit   none

      integer    size
      real*8     alpha , v1(*), v2(*) , w(*)

      integer    k
      real*8     a0

      a0 = 1.d0 - alpha
      do k = 1, size
        w(k) = a0*v1(k) + alpha*v2(k)
      end do ! k

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: bmatsvq

c     Purpose:    computes linearized strain matrix B

c     Input:      nel   : number of nodes in element
c     shp   :     shape functions and derivatives
c     phipr :     derivative w/r to S of centroid curve

c     Output:     bmat  : matrix B

c     Notes:      -
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine bmatsvq(nel, shp, phipr, bmat)

      implicit   none

      integer    nel
      real*8     shp(2,*), phipr(3), bmat(6,6,*)

      integer    a, i,j

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

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: checksvq

c     Purpose:    checks element for input errors in geometry

c     Input:      n   : node number
c     ndm :       spatial dimension
c     nel :       nodes in element
c     xl  :       reference nodal coordinates

c     Output:     only output messages

c     Notes :     only for 2 noded elements
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine checksvq(n, ndm, nel, xl)

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

        norm = 0.0d0
        do k = 1,ndm
          inc(k) = xl(k,2) - xl(k,1)
          norm   = norm + inc(k)*inc(k)
        end do ! k

        if ( norm .lt. 1d-5) then
          write (iow, 2000) n
          if (ior .lt. 0) then
            write(*, 2000) n
          endif
        endif

c       Second check, test for bad reference node/vector
        length = 1.d0/sqrt(norm)
        do k = 1,3
          t3(k) = inc(k)*length
        end do ! k

        if (lref.eq.1) then ! REVE NODE
          do k = 1,3
            t2(k) = refx(k) - xl(k,1)
          end do ! k
        elseif (lref.eq.2) then ! REVE VECTor
          do k = 1,3
            t2(k) = refx(k)
          end do ! k
        elseif (lref.eq.3) then ! REFE POLAr
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

c     Formats

 2000 format(//' Warning: Element ', i4, ' has zero length')
 2100 format(//' Warning: bad REFErence node/vector for element', i4)

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: commonssvq

c     Purpose: Obtains variables of commons /bmsvq1/,/bmsvq2/ from
c              element properties

c     Input:   d: material data

c     Output:  everything goes in common blocks

c     Notes:   At each integration point total rotation is stored. That
c     is rotation tensor to transforms inertial frame to section frame.
c     Observe there is also posibility to store incremental rotation
c     tensor that accounts for transformation between section frame
c     at time t=0 and current one.
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine commonssvq(d)

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

c     Select integration rule of stiffness (gint) & mass (mint) terms.
c     Use reduced integration for stiffness to avoid locking.

      if(nel.eq.2) then
        gint   = 1
        mint   = 2
      else
        gint   = 2
        mint   = 3
      endif
      if(nint(d(182)).gt.0) then      ! Nodal quadrature
        gint = mint
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

c     Set mct0, number of history terms for element. Per Gauss point:
c     In stiff: 4 (rotation quaternion)_n+3 spatial curvature omega_n
c     In mass : 4 (rotation quaternion)_n

      mct0 = 7*gint+4*mint

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
c     Subroutine: enersvq

c     Purpose:    computation of energy and momentum

c     Input:      nel: nodes in element
c                 ndf: degrees of freedom per node
c                 d  : feap material data
c                 xl : nodal reference positions
c                 ul : degrees of freedom vector
c                 isw: Switch value from element (should be 13)

c     Output:     In file P*.ene:

c                 epl(1-3) : linear momentum
c                 epl(4-6) : angular momentum
c                 epl(7)   : Kinetic energy
c                 epl(8)   : Potential energy
c                 epl(9)   : External energy
c                 epl(10)  : Total energy             (in pmacr2.f)
c                 epl(11)  : Angular momentum norm    (in pmacr2.f)

c     Notes:      for isw=13, independent of integrator, all variables
c                 (u,v,a) in ul() are at time tn+1
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine enersvq(nel , ndf , ndm , d , xl , ul, isw)

      implicit   none

      include   'bm3fd.h'         !mct0
      include   'cdata.h'         !nen
      include   'comblk.h'        !hr
      include   'erotas.h'        !rvel,
      include   'ddata.h'         !noi
      include   'hdata.h'         !nh1, nh2, nh
      include   'ptdat6.h'        !epl

      integer    nel         , ndf       , ndm , isw
      real*8     xl(ndm,*)   , ul(ndf,*) , d(*)

      logical    cons
      integer    nn1         , ip       , kk       , a    , nenv
      integer    nhs         , nns      , nh
      real*8     sg(2,5)     , shp(2,3) , xjac
      real*8     dx          , w
      real*8     xnode1(3,3) , omega1(3)
      real*8     phi1(3)     , phipr1(3)
      real*8     qrot1(4)
      real*8     rotmat1(3,3), strain(6)
      real*8     vn1(3)      , wn1(3)
      real*8     in1(3,3)    , pin1(3)   , tmp3(3)

c     flag that identifies energy-momentum integrator

      cons = noi .eq. 5

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

      nn1 = 0
      nh  = nint(d(15))    ! Number of history variables/point
      nhs = nint(d(149))   ! Number of history variables/section
      nns = mct0
      do ip = 1, gint

c       Recover data from history variables (rotation & curvature tn+1)

        do kk = 1,4
          qrot1(kk)  = hr(nh2+nn1-1+kk)
        end do ! kk
        do kk = 1,3
          omega1(kk) = hr(nh2+nn1+3+kk)
        end do ! kk

c       Shape functions and derivatives. Arclength element

        call shp1d(sg(1,ip), xl, shp, ndm, nel, xjac)
        dx = sg(2,ip)*xjac

c       Computation of position vectors and S-derivative at tn+1

        call positionsvq(xnode1, nel, shp, phi1, phipr1)

c       Computation of rotation matrix at time tn+1

        call quamat(qrot1 , rotmat1)

c       Compute convected strains.

        if(cons) then

c         Compute strain Gamma & Recover rotational strain from history

          do kk = 1,3
            strain(kk)   = rotmat1(1,kk)*phipr1(1)
     &                   + rotmat1(2,kk)*phipr1(2)
     &                   + rotmat1(3,kk)*phipr1(3)
            strain(3+kk) = omega1(kk)
          end do ! k
          strain(3) = strain(3) - 1.d0
         else

         endif

c        Compute stored energy

         call strenergysvq(d , hr(nh1+nns),hr(nh2+nns), nh,
     &                     strain, w, isw)

c        Accumulate stored enegy, integrated

         epl(8) = epl(8) + w*dx

         nn1    = nn1 + 7
         nns    = nns + nhs

      end do ! ip: Gauss loop

c     Gauss points and weights for full quadrature

      if(nint(d(182)).gt.0) then
        call int1dn(mint, sg)
      else
        call int1d(mint, sg)
      endif

      do ip = 1,mint

c       Shape functions and derivatives. Arclength element

        call shp1d(sg(1,ip), xl, shp, ndm, nel, xjac)
        dx = sg(2,ip)*xjac

c       Recover data from history variables (rotation_n+1)

        do kk = 1,4
          qrot1(kk)  = hr(nh2+nn1-1+kk)
        end do ! kk

c       Computation of position vectors and S-derivative at tn+1

        call positionsvq(xnode1, nel, shp, phi1, phipr1)

c       Rotation matrix at tn+1

        call quamat(qrot1 , rotmat1)

c       Linear and angular velocities at tn+1

        nenv = 3*nen
        do kk = 1,3
          vn1(kk) = 0.d0
          wn1(kk) = 0.d0
          do a = 1,nel
            vn1(kk) = vn1(kk) + shp(2,a) * ul(kk, nenv+a)
            wn1(kk) = wn1(kk) + shp(2,a) * rvel(kk,a,2)
          end do ! a
        end do ! kk

c       Angular momenta at time tn+1

        call momentsvq(rotmat1 , miner , wn1 , in1 , pin1)

c       Add up contributions for linear and angular momentum

        call vecp(phi1, vn1, tmp3)
        do kk = 1,3
          epl(kk)   = epl(kk)   + arho*vn1(kk)*dx
          epl(kk+3) = epl(kk+3) + arho*tmp3(kk)*dx + pin1(kk)*dx
        end do ! kk

c       Add up contributions for kinetic energy

        epl(7) = epl(7) + 0.5d0*dx*(arho*(vn1(1)*vn1(1)
     &                                  + vn1(2)*vn1(2)
     &                                  + vn1(3)*vn1(3))
     &                                  +(wn1(1)*pin1(1)
     &                                  + wn1(2)*pin1(2)
     &                                  + wn1(3)*pin1(3)))

c       Advance history variable pointers

        nn1 = nn1 + 4

      end do ! ip: Gauss loop

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: expmapsvq

c     Purpose:    compute exponential map of vector

c     Input:      v : (3) vector

c     Output:     r : (3,3) r = exp[v]
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine expmapsvq(v, r)

      implicit   none

      real*8     v(3) , r(3,3)
      real*8     q(4)

      call rotqua(v, q)
      call quamat(q, r)

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: forcesvq

c     Purpose:    compute distributed load

c     Input:      rotmatr: initial rotation tensor
c                 rotmat: rotation tensor
c                 eforce: body force, as input by user
c                 iforc : type of load.
c                   (1) : conservative load given in global axis
c                   (2) : conservative loald in rod local axis
c                   (3) : follower load in global axis
c                   (4) : follower load in rod local axis

c     Output:     bodyf : rotated body force
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine forcesvq(rotmatr, rotmat, eforce, iforc, bodyf)

      implicit   none

      integer    iforc
      real*8     rotmatr(3,3), rotmat(3,3), eforce(3), bodyf(3)

      integer    k
      real*8     tmp(3)

c     Conservative load in global axis

      if (iforc .eq. 1) then

        do k = 1,3
          bodyf(k) = eforce(k)
        end do ! k

c     Conservative load in local axis = Lambda_o*eforce

      else if (iforc .eq. 2) then

        do k = 1,3
          bodyf(k) = rotmatr(k,1)*eforce(1)
     &             + rotmatr(k,2)*eforce(2)
     &             + rotmatr(k,3)*eforce(3)
        end do ! k

c     Follower loads in global axis = Lambda*Lambda_o^T*eforce

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

c     Follower loads in local axis

      else if (iforc .eq. 4) then

        do k = 1,3
          bodyf(k) = rotmat(k,1)*eforce(1)
     &             + rotmat(k,2)*eforce(2)
     &             + rotmat(k,3)*eforce(3)
        end do ! k

      end if

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: hatsvq

c     Purpose:  Compute skew symmetric matrix associated to R3 vector

c     Input:    v  :(3)    vector

c     Output:   m  :(3,3)  matrix
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine hatsvq(v, m)

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
c  Subroutine: hmatsvq

c  Purpose:    Compute matrix h that appears when linearizing exp map

c  Input:      psi : rotation vector

c  Output:     h

c  Notes:      h = nn + sin(b)/b (eye-nn) + (1-cos(b))/n hatn,
c              b = norm(psi),
c              n = psi/b, nn= n otimes n
c              hmat = (tmat)^-1
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine hmatsvq(psi, h)

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
c     Subroutine: initializesvq

c     Purpose:    Initializes database at start of whole process

c     Input:      d   : feap material data
c                 xl  : reference nodal coordinates
c                 ndm : spatial dimension
c                 nel : nodes in element
c                 n   : node number

c     Output:     goes in common block

c     Notes:      The list of history variables initialized as follows:
c                 -- first all Gauss points for stiffness calculation
c                 -- then, all Gauss points for mass calculation
c                 For each G.p. for stiffness, following data stored
c                 -- a quaternion, with rotation
c                 -- spatial curvature, a (3) vector
c                 For each G.p. for mass:
c                 -- a quaternion, with rotatin

c                 Important: initialization routine commonssvq must be
c                            called before this one
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine initializesvq(d, xl, ndm, nel, n)

      implicit   none

      include   'bm3fd.h'
      include   'comblk.h'
      include   'erotas.h'
      include   'hdata.h'
      include   'refnd.h'
      include   'iofile.h'       !iow

      integer    ndm , nel , n
      real*8     d(*), xl(ndm, *)


      integer    ip, k, nn1, nn3
      integer    j1,k1,l1
      real*8     t1(3),t2(3),t3(3),length,sg(2,5),sgm(2,5)
      real*8     xlbda0(3,3),shp(2,3),xjac,aa

      save

c     External node, and coordinates
c     In material definition, one must give one of following
c     REFE, NODE , x1, x2, x3 ... and coordinates of node  t1, or
c     REFE, VECT,  v1, v2, v3 ... and components of vector t1
c     REFE, POLA,

      lref    = nint(d(96))
      refx(1) = d(97)
      refx(2) = d(98)
      refx(3) = d(99)

c     Checking of a straight element

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

        if (ior .lt. 0)
     &       write(*,*) 'Warning: bad REFErence  for element' ,n

      else

        length = 1.d0/length
        do j1    = 1,3
          t1(j1) = t1(j1)*length
         end do ! j1
         call vecp(t3,t1,t2)

      end if

c     Form rotation matrix and extract quaternion

      do j1 = 1,3
        xlbda0(j1,1) = t1(j1)
        xlbda0(j1,2) = t2(j1)
        xlbda0(j1,3) = t3(j1)
      end do ! j1

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

 101  if(nint(d(182)).gt.0) then
        call int1dn(gint, sg)
      else
        call int1d(gint, sg)
      endif

      do j1 = 1,gint
        call shp1d(sg(1,j1), xl, shp, ndm, nel, xjac)
        do k1 = 1,3
          t3(k1) = 0.0d0
        end do                 ! k1

c       Third vector has direction of tangent to line of centroids

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
        call int1d (mint,sgm)
      else
        call int1dl(mint,sgm)
      end if
      do j1 = 1,mint
        call shp1d(sgm(1,j1),xl,shp,ndm,nel,xjac)
        do k1 = 1,3
          t3(k1) = 0.0d0
        end do ! k1

c       Third vector has direction of tangent to line of centroids

        do l1 = 1,nel
          do k1 = 1,3
            t3(k1) = t3(k1) + shp(1,l1)*xl(k1,l1)
          end do !k1
        end do !l1
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

c     Rotation matrix is initialized for all stiffness GP of element

 200  nn1 = 0
      do ip = 1,gint

        do k = 1,4
          hr(nh1+nn1-1+k) = rotqu0(k, ip)
          hr(nh2+nn1-1+k) = rotqu0(k, ip)
        end do ! k

c       Set initial curvature to 0
        do k = 1,3
          hr(nh1+nn1+3+k) = 0.d0
          hr(nh2+nn1+3+k) = 0.d0
        end do ! k

        nn1 = nn1 + 7
      end do ! ip

c     Rotation matrix is initialized for all mass GP of element

      do ip = 1,mint

        do k = 1,4
          hr(nh1+nn1-1+k) = rotqu0(k, gint+ip)
          hr(nh1+nn1-1+k) = rotqu0(k, gint+ip)
        end do ! k

        nn1 = nn1 + 4
      end do ! ip

c     This rotation matrices are stored in nh3, they never change
c     because we will need to know initial rotation

      nn3 = 0
      do ip = 1,gint+mint

        do k = 1,4
          hr(nh3+nn3-1+k) = rotqu0(k,ip)
        end do ! k
        nn3 = nn3 + 4

      end do ! ip

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: kgeosvq

c     Purpose:    compute geometric part of tangent stiffness

c     Input:      nel : nodes on element
c                 shp : shape functions and derivatives
c               phipr : derivative w/r to S of centroid curve
c                thts : rotation vector
c              sstress: spatial stresses [n m]

c     Output:     kgeo : tangent
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine kgeosvq(nel, shp, phipr, sstress, kgeo)

      implicit   none

      integer    nel
      real*8     shp(2,*), phipr(3), sstress(6)
      real*8     kgeo(18,18)

      integer    a, b, aa, bb , i, j
      real*8     ndphi
      real*8     hatn(3,3) , hatm(3,3)

      call hatsvq( sstress(1), hatn )
      call hatsvq( sstress(4), hatm )
      ndphi = sstress(1)*phipr(1)
     &      + sstress(2)*phipr(2)
     &      + sstress(3)*phipr(3)

      aa = 0
      do a = 1,nel
        bb = 0
        do b = 1,nel

c         Block 1-1

          do j = 1,3
            do i = 1,3
              kgeo(aa+i,bb+j) =  0.0d0
             end do ! i
          end do ! j

c         Block 1-2

          do j = 1,3
            do i = 1,3
              kgeo(aa+i,bb+3+j) =  -shp(1,a)*shp(2,b)*hatn(i,j)
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
              kgeo(aa+3+i,bb+3+j) = (shp(2,a)*sstress(i)*phipr(j)
     &                            -  shp(1,a)*hatm(i,j))*shp(2,b)
            end do ! i
            kgeo(aa+3+j,bb+3+j) = kgeo(aa+3+j,bb+3+j)
     &                          - shp(2,a)*shp(2,b)*ndphi
          end do ! j
          bb = bb+6
        end do ! bb
        aa = aa+6
      end do ! aa

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: kmatsvq

c     Purpose:    compute material part of tangent stiffness

c     Input:
c                  nel : nodes on element
c                  bmat: linearized strain matrix at equilibrium point
c                  scc : spatial material elasticity

c     Output:     kmat : tangent
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine kmatsvq(nel , bmat , scc, kmat)

      implicit   none

      integer    nel
      real*8     bmat(6,*), scc(6,6), kmat(18,18)

      integer    j,k,jj,kk
      real*8     btc(6,18)

c     Multiply Bt C and store its transpose

      jj = 0
      do j = 1,nel
        do k = 1,6
          btc(k,jj+1) = bmat(1,jj+1)*scc(1,k)
          btc(k,jj+2) = bmat(2,jj+2)*scc(2,k)
          btc(k,jj+3) = bmat(3,jj+3)*scc(3,k)
          btc(k,jj+4) = bmat(1,jj+4)*scc(1,k) + bmat(2,jj+4)*scc(2,k)
     &                + bmat(3,jj+4)*scc(3,k) + bmat(4,jj+4)*scc(4,k)
     &                + bmat(5,jj+4)*scc(5,k) + bmat(6,jj+4)*scc(6,k)
          btc(k,jj+5) = bmat(1,jj+5)*scc(1,k) + bmat(2,jj+5)*scc(2,k)
     &                + bmat(3,jj+5)*scc(3,k) + bmat(4,jj+5)*scc(4,k)
     &                + bmat(5,jj+5)*scc(5,k) + bmat(6,jj+5)*scc(6,k)
          btc(k,jj+6) = bmat(1,jj+6)*scc(1,k) + bmat(2,jj+6)*scc(2,k)
     &                + bmat(3,jj+6)*scc(3,k) + bmat(4,jj+6)*scc(4,k)
     &                + bmat(5,jj+6)*scc(5,k) + bmat(6,jj+6)*scc(6,k)
        end do ! k
        jj = jj + 6
      end do ! j

      do j = 1, 6*nel
        kk = 0
        do k = 1,nel
          kmat(j,kk+1) = btc(1,j)*bmat(1,kk+1)
          kmat(j,kk+2) = btc(2,j)*bmat(2,kk+2)
          kmat(j,kk+3) = btc(3,j)*bmat(3,kk+3)
          kmat(j,kk+4) = btc(1,j)*bmat(1,kk+4) + btc(2,j)*bmat(2,kk+4)
     &                 + btc(3,j)*bmat(3,kk+4) + btc(4,j)*bmat(4,kk+4)
     &                 + btc(5,j)*bmat(5,kk+4) + btc(6,j)*bmat(6,kk+4)
          kmat(j,kk+5) = btc(1,j)*bmat(1,kk+5) + btc(2,j)*bmat(2,kk+5)
     &                 + btc(3,j)*bmat(3,kk+5) + btc(4,j)*bmat(4,kk+5)
     &                 + btc(5,j)*bmat(5,kk+5) + btc(6,j)*bmat(6,kk+5)
          kmat(j,kk+6) = btc(1,j)*bmat(1,kk+6) + btc(2,j)*bmat(2,kk+6)
     &                 + btc(3,j)*bmat(3,kk+6) + btc(4,j)*bmat(4,kk+6)
     &                 + btc(5,j)*bmat(5,kk+6) + btc(6,j)*bmat(6,kk+6)
          kk = kk + 6
        end do ! k
      end do ! j

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: masssvq

c     Purpose:    computation of mass matrix

c     Input:      xl : nodal reference positions
c                 nel: nodes in element
c                 ndm: spatial dimension
c                 nst: dimension of mass matrix (ndm*nel)

c     Output:     s: mass matrix
c                 p: lumped mass matrix

c     Notes:      mass matrix is not consistant for this model
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine masssvq(xl , nel , ndf , ndm , nst , s , p)

      implicit   none

      include   'bm3fd.h'
      include   'comblk.h'
      include   'cdata.h'
      include   'ddata.h'         !noi, theta
      include   'eltran.h'
      include   'erotas.h'
      include   'hdata.h'
      include   'tdata.h'

      integer    ndf        , ndm        , nst       , nel
      real*8     xl(ndm,*)  , s(nst,nst) , p(nst)

      logical    cons
      integer    a, b, i, k1,l1, kk, ip
      integer    nn1

      real*8     dx         , xjac
      real*8     sg(2,5)    , shp(2,3)
      real*8     qrot0(4)   , qrot(4)     , qrot1(4)
      real*8     thts(3)    , tht1(3)
      real*8     rotmat1(3,3)
      real*8     wn1(3)     , pin1(3)     , in1(3,3)
      real*8     hpin1(3,3)
      real*8     tmatrix(3,3)
      real*8     tmp        , tmp2(3,3)

      cons = noi .eq. 5

c     Computation of Gauss points and weights for stiffness integration

      call int1d (mint , sg)

c     Loop on each Gauss point

      nn1   = 7*gint
      do ip = 1,mint

c       Recover data from history variables (rotation_n)

        do kk = 1,4
          qrot0(kk)  = hr(nh1+nn1-1+kk)
        end do ! kk

c       Shape functions and derivatives. Arclength element

        call shp1d(sg(1,ip), xl, shp, ndm, nel, xjac)
        dx = sg(2,ip)*xjac

c       Spatial rotation vector from Lambda_n to Lambda_n+alpha

        do kk = 1,3
          thts(kk) = 0.0d0
          tht1(kk) = 0.0d0  ! ERROR Missing ?
          do a = 1,nel
            thts(kk) = thts(kk) + theta(3)*rots(kk,a,2)*shp(2,a)
            tht1(kk) = tht1(kk) +          rots(kk,a,2)*shp(2,a)
          end do ! a
        end do ! kk

c       Computation of rotation_n+1, rotation_n+a at Gauss point
c       (quaternion for

        call updrotationsvq(qrot0 , tht1 , qrot1)
        call updrotationsvq(qrot0 , thts , qrot )

c       Computation of rotation tensor (matrices) at time tn and tn+1

        call quamat(qrot1 , rotmat1)

c       Angular velocities

        do i = 1,3
          wn1(i) = 0.d0
          do a = 1,nel
            wn1(i) = wn1(i) + shp(2,a) * rvel(i,a,2)
          end do ! a
        end do ! i

c       Angular momenta at time  tn+1

        call momentsvq(rotmat1 , miner , wn1 , in1 , pin1)

c       Calculations for rotational block

        call hatsvq(pin1 , hpin1)
        call tmatsvq(thts, tmatrix)

c       multiply spatial inertia (symmetric) by tmatrix^T

        do b = 1,3
          do a = 1,3
            tmp2(a,b) = in1(1,a)*tmatrix(b,1)
     &                + in1(2,a)*tmatrix(b,2)
     &                + in1(3,a)*tmatrix(b,3)
          end do ! a
        end do ! b


        l1 = 0
        do a = 1,nel

          k1 = 0
          do b = 1,nel

c           Translation contribution

            tmp          = shp(2,a)*shp(2,b)*dx*ctan(3)*arho

            s(l1+1,k1+1) = s(l1+1,k1+1) + tmp
            s(l1+2,k1+2) = s(l1+2,k1+2) + tmp
            s(l1+3,k1+3) = s(l1+3,k1+3) + tmp

c           Rotation contribution

c           tmp          = shp(2,a)*shp(2,b)*dx*ctan(3)

c           do i = 1,3
c             do kk = 1,3
c               s(l1+i+3,k1+kk+3) = s(l1+i+3,k1+kk+3) + 1.0d0
c     &                           + tmp2(i,kk)*ctan(3)
c     &                           + hpin1(kk,i)/(dt*ctan(3))
c             end do ! j
c           end do ! i

            k1 = k1 + ndf
          end do ! b
          l1 = l1 + ndf
        end do ! a

c       Lumped mass matrix, just add up columns

        do i = 1,nst
           do a = 1,nst
              p(a) = p(a) + s(i,a)
           end do ! a
        end do ! i

c       advance history variable pointers

        nn1 = nn1 + 4

      end do ! ip: Gauss loop

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: modelsvq

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
      subroutine modelsvq(d , hn, h1, nh, strain, stress, cc, isw)

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
c     Subroutine: momentsvq

c     Purpose:    compute spatial inertia and angular momentum

c     Input:      rot  : (3,3) rotation tensor
c                 miner: (3,3) body inertia
c                 w    : (3)   spatial angular velocity

c     Output:     in   : (3,3) spatial inertial
c                 pi   : (3)   spatial section momentum

c     Notes:      pi = rot * miner * rot^T * w
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine momentsvq(rot , miner, w, in , pi)

      implicit   none

      real*8     rot(3,3) , miner(3,3) , w(3) , in(3,3) , pi(3)

      integer    i , j
      real*8     tmp(3,3)

c     Compute spatial inertia by pushing forward material inertia
c     i = Lambda * I * Lambda^T

      do j = 1,3
        do i = 1,3
          tmp(i,j) = miner(i,1)*rot(j,1)
     &             + miner(i,2)*rot(j,2)
     &             + miner(i,3)*rot(j,3)
        end do ! i
      end do ! j

      do j = 1,3
        do i = 1,3
          in(i,j) = rot(i,1)*tmp(1,j)
     &            + rot(i,2)*tmp(2,j)
     &            + rot(i,3)*tmp(3,j)
        end do ! i
      end do ! j

c     Spatial momentum

      do i = 1,3
        pi(i) = in(i,1)*w(1) + in(i,2)*w(2) + in(i,3)*w(3)
      end do ! i

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: positionsvq

c     Purpose:    Computes current position of one point on line of
c                 centroids (Phi) and S-derivative (Phi,s)

c     Input:      xnode : (3,nel) nodal position vectors
c                 nel   : number of nodes in element
c                 shp   : (2,nel) derivative and shape functions

c     Output:     phi   : (3) position vector at Gauss point
c                phipr  : (3) derivative of phi w/r to reference
c                         arclength
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine positionsvq(xnode , nel , shp , phi , phipr)

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
c     Subroutine: pushfsvq

c     Purpose:    compute algorithmic spatial stresses and material
c                 elasticities by pushing forward their corresponding
c                 convective tensors

c     Input:      rotmat0: (3,3) rotation matrix at time tn
c                 stress : (6)   convected stresses [N , M]
c                 cc     : (6,6) convected elasticities

c     Output:     sstress: (6)   spatial stresses [n , m]
c                 scc    : (6,6) spatial elasticities
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine pushfsvq(rotmat, stress, cc, sstress , scc)

      implicit   none

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
c     Subroutine: statsvq

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
      subroutine statsvq(d, xl, ul, s, p, ndf ,ndm, nst, isw, xnode)

      implicit   none

      include   'bm3fd.h'
      include   'cdata.h'
      include   'comblk.h'
      include   'counts.h'
      include   'ddata.h'        !noi
      include   'debugs.h'
      include   'eldata.h'
      include   'elplot.h'       !tt
      include   'eltran.h'
      include   'erotas.h'       ! xln, rots
      include   'fdata.h'        ! fl
      include   'hdata.h'        ! nh1, nh2
      include   'iofile.h'       ! iow, ior
      include   'tdata.h'

      integer    ndf      , ndm       , nst      , isw
      real*8     d(*)     , xl(ndm,*) , ul(ndf,*), s(nst,nst), p(nst)
      real*8     xnode(3,*)


      integer    a, b, k1,l1, kk, ip
      integer    nn1, nn3, nns, nh, nhs

      real*8     dx        , xjac
      real*8     sg(2,5)   , shp(2,3)
      real*8     phi(3)    , phipr(3)    , pp(6,2)
      real*8     qrotr(4)  , qrot0(4)    , qrot(4)
      real*8     omega0(3) , omega(3)
      real*8     thts(3)   , thtp(3)
      real*8     rotmat(3,3)
      real*8     rotmatr(3,3)
      real*8     bodyf(3)
      real*8     strain(6)
      real*8     stress(6) , sstress(6)
      real*8     cc(6,6)   , scc(6,6)
      real*8     bmat(6,6,3)
      real*8     kmat(18,18), kgeo(18,18)
      real*8     tmp

      logical    dynamics, cons

      save

c     flags, dynamics->true when running dynamics,
c                cons->true when trans,cons

      dynamics = fl(9)
      cons     = noi .eq. 5

c     Computation of Gauss points and weights for stiffness integration

      call int1d (gint , sg)

c     Loop on each Gauss point

      nn1 = 0
      nn3 = 0
      nns = mct0
      nh  = nint(d(15))  ! Number of history variables/quad. point
      nhs = nint(d(149)) ! Number of history variables/section
      do ip = 1,gint

c       Recover data from history variables (rotation_n, rotation_ref)

        do k1 = 1,4
          qrot0(k1)  = hr(nh2+nn1-1+k1)
          qrotr(k1)  = hr(nh2+nn3-1+k1)
        end do ! k1

c       Recover convected rotational strain at t_n

        do k1 = 1,3
          omega0(k1) = hr(nh2+nn1+3+k1)
        end do ! k1

c       Shape functions and derivatives. Arclength element

        call shp1d(sg(1,ip), xl, shp, ndm, nel, xjac)
        dx = sg(2,ip)*xjac

c       Computation of position vectors & derivative with respect to S

        call positionsvq(xnode , nel, shp, phi , phipr)

c       Compute spatial rotation vector that takes Lambda_n to
c       Lambda_n+1.  Also its derivative, for computation of curvature

        do kk = 1,3
          thts(kk) = 0.0d0
          thtp(kk) = 0.0d0
          do a = 1,nel
             thts(kk) = thts(kk) + rots(kk,a,1)*shp(2,a)
             thtp(kk) = thtp(kk) + rots(kk,a,1)*shp(1,a)
          end do ! a
        end do ! kk

c       Computation of rotation_n+1 at Gauss point (quaternion format)

        call updrotationsvq(qrot0, thts, qrot)

c       Computation of rotation tensor (matrices) at time t=0, tn
c       and tn+1

        call quamat(qrotr , rotmatr)
        call quamat(qrot  , rotmat)

c       convected strain vector omega at tn+1

        call updomegasvq(omega0, thts, thtp, rotmat, omega)

c       strain tensor and linearized strain matrix

        call strainsvq(phipr, rotmat , omega, strain)
        call bmatsvq(nel, shp, phipr, bmat)

c       Compute convected stresses & material tangent: push to spatial

        call modelsvq(d, hr(nh1+nns), hr(nh2+nns), nh,
     &                strain, stress, cc, isw)
        call pushfsvq(rotmat, stress, cc, sstress, scc)

c       Compute contribution from distributed forces (conservative or
c       follower) and contribution to tangent

        call forcesvq(rotmatr, rotmat, eforce, iforc, bodyf)

c       Static contribution to residual

        k1 = 0
        do a = 1,nel

c         stress divergence terms

          do kk = 1,6
            p(k1+kk) = p(k1+kk) - (bmat(1,kk,a)*sstress(1)
     &                           + bmat(2,kk,a)*sstress(2)
     &                           + bmat(3,kk,a)*sstress(3)
     &                           + bmat(4,kk,a)*sstress(4)
     &                           + bmat(5,kk,a)*sstress(5)
     &                           + bmat(6,kk,a)*sstress(6))*dx
          end do ! kk

c         body forces (no body couples) contribution

          do kk = 1,3
            p(k1+kk) = p(k1+kk) + shp(2,a)*bodyf(kk)*dx
          end do ! kk
          k1  = k1 + ndf
        end do ! a

c       Computation of tangent

        if (isw .eq. 3) then

          dx = dx * ctan(1)

c         Material part of tangent
          call kmatsvq(nel, bmat, scc, kmat)

c         Geometric part
          call kgeosvq(nel, shp, phipr, sstress, kgeo)

c         Assemble both contributions in element tangent
          do l1 = 1,nel*6
            do k1 = 1,nel*6
              s(is(l1),is(k1)) = s(is(l1),is(k1))
     &                         + (kgeo(l1,k1)+kmat(l1,k1))*dx
            end do ! k1
          end do ! i1

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
          write(iow,2002) strain, stress

          if(ior.lt.0) then
            write(*,2001) n,ip,phi,rotmat
            write(*,2002) strain, stress
          endif

        end if

c       Store data to history variables (rotation_n+1 & rot.strain_n+1)

        do k1 = 1,4
          hr(nh2+nn1-1+k1) = qrot(k1)
        end do ! k1

        do k1 = 1,3
          hr(nh2+nn1+3+k1) = omega(k1)
        end do ! k1

        nn1 = nn1 + 7
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

c     Project stresses to nodes for plots

      if(isw .eq. 8) then
        do kk = 1,6
          pp(kk,1) =  p(kk)
          pp(kk,2) = -p(kk+ndf)
        end do ! kk
        call frcn3d(pp,p,s)
      endif

c     Formats

 2001 format(5x,'Element No.',i7,': Point No.',i2/
     &     5x,'Current Coordinates'/5x,1p,3e15.5/
     &     5x,'Current Intrinsic Frame'/5x,3(1p,3e15.5/5x))

 2002 format(5x,'Material Strains'/
     &     5x,'    Gamma : ',1p,3e15.5/
     &     5x,'    Omega : ',1p,3e15.5/
     &     5x,'Material Stresses'/
     &     5x,'    N     : ',1p,3e15.5/
     &     5x,'    M     : ',1p,3e15.5/)

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: strainsvq

c     Purpose:    Computes convected strain tensors

c     Input:      phipr    : derivative w/r to S of centroid curve
c                 rotmat  : (3,3) rotation matrix
c                 omg     : (3) convected rotational strain

c     Output:     strain  : (6) convected strain [Gamma Omega]

c     Notes:      Only valid for straight elements in reference
c                 configuration. Only modification needed is
c                 to substract ref. strains from current ones.
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine strainsvq(phipr, rotmat, omg , strain)

      implicit   none

      integer    i
      real*8     phipr(3), rotmat(3,3)
      real*8     omg(3)
      real*8     strain(6)

      save

c     Computation of strain tensor

      do i = 1,3
         strain(i)   = rotmat(1,i)*phipr(1)
     &               + rotmat(2,i)*phipr(2)
     &               + rotmat(3,i)*phipr(3)
         strain(i+3) = omg(i)
      end do ! i
      strain(3) = strain(3) - 1.d0

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: strenergysvq

c     Purpose:    Computes stored energy function

c     Input:      d       : feap's material data
c                 strain  : convected strains
c                   1-3   : strain Gamma
c                   4-6   : strain Omega
c                 isw     : Switch option

c     Output:     w       : stored energy function

c     Notes:      Only valid for quadratic potentials, using model
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine strenergysvq(d , hn, h1, nh, strain , w, isw)

      implicit   none

      integer    nh   , isw
      real*8     d(*) , hn(*), h1(*), strain(6)
      real*8     w

      real*8     stress(6) , cc(6,6)

      call modelsvq(d, hn, h1, nh, strain, stress, cc, isw)
      w = 0.5d0*(strain(1)*stress(1) + strain(2)*stress(2)
     &         + strain(3)*stress(3) + strain(4)*stress(4)
     &         + strain(5)*stress(5) + strain(6)*stress(6))

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: tangsvq

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
      subroutine tangsvq(d, ul, xl, s, p, ndf, ndm, nst, isw, symflg)

      implicit   none

      include   'bm3fd.h'
      include   'cdata.h'
      include   'ddata.h'
      include   'eldata.h'
      include   'eltran.h'
      include   'erotas.h'       !xln
      include   'fdata.h'        !fl
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

      call statsvq(d, xl, ul, s, p, ndf, ndm , nst, isw, xnode)

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
c     Subroutine: tmatsvq

c     Purpose:    compute matrix T in linearization of rotation vector

c     Input:      thts : rotation vector

c     Output:     tmatrix

c     Notes:      Different if updates are done with exponential map
c                 or Cayley transform
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine tmatsvq(thts, tmatrix)

      implicit   none

      include   'erotas.h'        !rotyp

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
c     Subroutine: updomegasvq

c     Purpose:    Update spatial curvature omega

c     Input:
c                 omega0: (3) spatial curvature vector omega at tn
c                 tht   : (3) spatial rotation vector
c                 thtp  : (3) derivative of spatial rotation vector
c                 rotmat: (3x3)rotation matrix at equilibrium.

c     Output:     omega1: (3) updated spatial curvature vector

c     Notes:      Implementation equivalent to Simo & Vu Quoc[1986]
c                 given in appendix B. They update spatial curvature,
c                 convected strain is updated here.
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine updomegasvq(omega0 , tht , thtp , rotmat, omega1)

      implicit   none

      real*8     omega0(3),  tht(3), thtp(3), rotmat(3,3) , omega1(3)

      integer    i
      real*8     hmat(3,3), tmp(3)

      call hmatsvq(tht, hmat)

      do i=1,3
        tmp(i) = hmat(i,1)*thtp(1)
     &         + hmat(i,2)*thtp(2)
     &         + hmat(i,3)*thtp(3)
      end do ! i

      do i=1,3
        omega1(i) = omega0(i) + rotmat(1,i)*tmp(1)
     &                        + rotmat(2,i)*tmp(2)
     &                        + rotmat(3,i)*tmp(3)
      end do ! i

      end

c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
c     Subroutine: updrotationsvq

c     Purpose: Update rotation matrix given spatial incremental
c              rotation vector,
c              using exponential map -> rot2 = exp[tht] rot1
c              or Cayley transform   -> rot2 = cay[tht] rot1

c     Input:      cons  : true if conserving method, so use cayley
c                 rot   : (4) rotation tensor in quaternion format
c                 tht   : (3) rotation vector

c     Output:     rot2  : (4) updated rotation tensor, quaternion format
c-----[--+----+----+----+----+----+----+----+----+----+----+----+----+-]
      subroutine updrotationsvq(qrot0 , tht , qrot1)

      implicit   none

      real*8     qrot0(4), tht(3), qrot1(4), qtht(4)

c     Computation of quaternions associated to tht

      call rotqua(tht, qtht)

c     Computation of quaternions of new rotation matrix

      call quamul(qtht, qrot0, qrot1)
      call quanrm(qrot1)

      end
