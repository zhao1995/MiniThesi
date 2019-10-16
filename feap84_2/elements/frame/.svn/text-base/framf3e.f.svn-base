c$Id:$
      subroutine framf3e(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove mprint from debug print                   04/10/2008
c-----[--.----+----.----+----.-----------------------------------------]
c     3-D Finite Deformation Frame (INCLUDING UPDATES FOR SO(3) AND S2)

c       3-D Geometrically-exact beam
c       Linear elastic constitutive equations
c       Finite displacements and strains
c       Nonlinear dynamic analysis by Newton-Raphson method
c       Uniformly distributed non-follower loads and no external torques
c       2 nodes (linear) or 3 nodes (quadratic) straight element

c       Total updates of variables from t=n to t=n+alpha
c       Independent variables: Incremental translational displacements
c       and rotation vectors from t=n to t=n+alpha
c       Interpolation of linear and angular momenta, spatial angular
c       velocities and rotation vector

c       J. C. Simo,  M. Doblare    -     July 1990

c       Modified to Energy-Momentum method with Cayley transform
c       N. Tarnow                  -     January 1993

c       Modified to standard FEAP element
c       R. Taylor                  -     January 1997
c-----[--.----+----.----+----.-----------------------------------------]
c   Partition nh1

c    xlbds(4,gint) - Quaternions associated with rotation matrices at
c                    Gauss points for stiffness at t=n
c    xlbdd(4,mint) - Quaternions associated with rotation matrices at
c                    Gauss points for mass at t=n
c    omega(3,gint) - Spatial curvatures at Gauss points for stiffness
c                    at t=n

c   Partition nh2

c    thts(3,gint)  - Spatial rotation vectors at Gauss points of stiff
c                    between n and n+alpha
c    thtpm(3,gint) - Derivative with respect to S of spatial rotation
c                    vectors at Gauss points for stiffness between n and
c                    n+alpha
c    thtd(3,gint)  - Spatial rotation vectors at Gauss points for mass
c                    between n and n+alpha

c   /bm3f1/  - One of internal Commons of routine FRM3FD
c    gint          - Number of Gauss points for stiffness
c    lobatt        - Flag that indicates if a Gauss-Legendre quadrature
c                    must be used for mass matrix (0) or a Gauss-Lobatto
c                    one (1)
c    mint          - Number of integration points for mass
c    isym          - Flag that indicates if tangent stiffness must be
c                    symmetrized (1) or not (0)
c    iforc         - Flag indicating type of external force:
c                    1 - nonfollower forces in global coordinates
c                    2 - nonfollower forces in local coordinates
c                    3 - follower forces along directors
c                    4 - follower forces along line of centroid
c    mct0          - Maximum of sizes of partition of array "h"
c                    (nh1 and nh2) used by this routine

c /bm3f2/  - One of internal Commons of routine FRM3FD
c    arho          - The section area weighted by reference density
c    miner(3,3)    - Material inertia dyadic
c    mhook(6,6)    - Material hookean tensor
c    eforce(3)     - Distributed loads along element at t=n+alfa
c    rotqu0(4,13)  - Initial rotation quaternion at Gauss points for
c                    mass and stiffness and nodes
c-----[--.----+----.----+----.-----------------------------------------]
c CALL : INIT3F, INPD3F, TSRF3F, OUTP3F, INDB3F
c        UPDB3F, TSDB3F, PRCM3F, PRSP3F
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3cn.h'
      include  'bm3fd.h'
      include  'comblk.h'
      include  'cdata.h'
      include  'ddata.h'
      include  'eldata.h'
      include  'erotas.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'ptdat6.h'
      include  'pview.h'
      include  'tdata.h'

      integer   ndf,ndm,nst,isw,ix(3)
      integer   i, k2,k3,k5,k6,nh
      real*8    xlbd1(4),d(*),ul(ndf,*),xl(ndm,nel),s(nst,nst),p(nst)

      save

      if(isw.eq.0) then
        if(ior.lt.0) then
        write(*,*) '   FRM3FD: 3-d Frame Finite Defm (2/3 node)'
        endif
        return

c     Some previous checkings

      elseif(isw.gt.2) then

c       Initial rotation matrices for element

        call inir3f(d,xl,ndm,nel)

c       Pointers for variables in partition of array 'h'

        k2 = 4*gint
        k3 = k2+4*mint

        k5 = 3*gint
        k6 = k5+3*gint

        nh  = int(d(15))   ! Number of history terms/stress point

      endif

c     Input element parameters

      if(isw.eq.1) then

c       Initialization of common /bm3f1/, /bm3f2/ and computation of mct

        call init3f(d)
        nh1 = nh1 + mct0

c       Set rotational update type

        rotyp = -4

c     Computation of tangent stiffness and residual vector

      elseif(isw.eq.3 .or. isw.eq.6) then

        call init3f(d)
        call tsrf3f(d,ul,xl,s,p,ndf,ndm,nst,hr(nh1),hr(nh1+k3),
     &              hr(nh2),hr(nh2+k5),hr(nh1+mct0),hr(nh2+mct0),
     &              nh,isw)

c     Output of variables

      elseif(isw.eq.4 .or.isw.eq.8) then

        call init3f(d)
        call outp3f(d,ul,xl,s,p,ndf,ndm,nst,hr(nh1),hr(nh1+k3),
     &              hr(nh2),hr(nh2+k5),hr(nh1+mct0),hr(nh2+mct0),
     &              nh,isw)

c     Compute element mass matrices

      elseif(isw.eq.5) then

        call init3f(d)
        call mass3f(xl,s,p,ndf,ndm,nst)

c     Computation of projected surface stresses

      elseif(isw.eq.20) then

        call init3f(d)
        do i = 1,nel
          mxl(i) = ix(i)
          xll(1,i) = xl(1,i)
          xll(2,i) = xl(2,i)
          xll(3,i) = xl(3,i)
        end do ! i
        if(cs.gt.0.0d0) then
          do i = 1,nel
            xll(1,i) = xll(1,i) + cs*ul(1,i)
            xll(2,i) = xll(2,i) + cs*ul(2,i)
            xll(3,i) = xll(3,i) + cs*ul(3,i)
          end do ! i
          call uprm3f(hr(nh1),hr(nh2),xlbd1)
          call quamat(xlbd1,rot1)
        endif
        mel    = nel
        mxl(4) = ma

c     Update of database just before start of a new time step

      elseif(isw.eq.12) then

        call init3f(d)
        call tsdb3f(hr(nh1),hr(nh1+k2),hr(nh1+k3),
     &              hr(nh2),hr(nh2+k5),hr(nh2+k6))

c     Computes and prints angular momentum and kinetic, potential
c     and total energies

      elseif(isw.eq.13) then

        if(theta(1).ne.0.d0.and.ttim.ne.0.d0) then
          call init3f(d)
          call agen3f(d,ul,xl,ndf,ndm,hr(nh1),hr(nh1+k3),
     &                hr(nh2),hr(nh2+k5),
     &                hr(nh1+mct0),hr(nh2+mct0),nh,isw)
        endif

c     Initialize data base values

      elseif(isw.eq.14) then

        call indb3f(hr(nh1),hr(nh1+k2),hr(nh1+k3),
     &              hr(nh2),hr(nh2+k5),hr(nh2+k6))

c     Other values of "isw" are not used in this routine

      endif

      end

      subroutine tsrf3f(d,ul,xl,s,p,ndf,ndm,nst,xlbds,omega,thts,thtpm,
     &            hn,h1,nh,isw)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Computes tangent stiffness and residual forces
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i o   xlbds(4,gint)  : Quaternions associated with rotation matrices
c                        at Gauss points for stiffness at t=n
c i o   omega(3,gint)  : Curvatures at Gauss points for stiffness at t=n
c i o   thts(3,gint)   : Current spatial rotation vectors at Gauss point
c                        for stiffness between t=n and t=n+alpha
c i o   thtpm(3,gint)  : Current derivative with respect to S of spatial
c                        rotation vectors at Gauss points for stiffness
c                        between t=n and t=n+alpha
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: FRM3FD
c CALL     : INDB3F, UPDB3F, STAT3F, DYNA3F, PRSP3F
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'
      include  'cdata.h'
      include  'ddata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'part0.h'
      include  'rdata.h'
      include  'tdata.h'

      integer   ndf,ndm,nst,nh,isw, j1,k1
      real*8    d(*),ul(ndf,*),xl(ndm,nel),s(nst,nst),p(nst),hn(*),h1(*)
      real*8    omega(3,gint),thts(3,gint),thtpm(3,gint),xlbds(4,gint)

      save

c     Update database

      call updb3f(d,xl,ndm,thts,thtpm)

c     Compute static  tangent stiffness and residual vector

      call stat3f(d,ul,xl,s,p,ndf,ndm,nst,xlbds,omega,thts,thtpm,
     &            hn,h1,nh,isw)

c     Compute dynamic tangent stiffness and residual vector

      if(ctan(3).ne.0.0d0 .and. ((ndfo(1).gt.0 .and. dt.gt.0.0d0) .or.
     &                                                  shflg) ) then
        call dyna3f(ul,xl,s,p,ndf,ndm,nst)
      end if

c     Symmetrization of stiffness matrix

      if(isym.eq.1) then

        do j1 = 1,nst
          do k1 = j1,nst
            s(j1,k1) = (s(j1,k1) + s(k1,j1))*0.5d0
          end do ! k1
        end do ! j1
        do j1 = 1,nst
          do k1 = 1,j1-1
            s(j1,k1) = s(k1,j1)
          end do ! k1
        end do ! j1
      end if

      end

      subroutine stat3f(d,ul,xl,s,p,ndf,ndm,nst,xlbds,omega,thts,thtpm,
     &                  hn,h1,nh,isw)

c     Computes static part of tangent stiffness and residual forces
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i     xlbds(4,gint)  : Quaternions associated with rotation matrices
c                        at Gauss points for stiffness at t=n
c i     omega(3,gint)  : Curvatures at Gauss points for stiffness at t=n
c i o   thts(3,gint)   : Current spatial rotation vectors at Gauss point
c                        for stiffness between t=n and t=n+alpha
c i o   thtpm(3,gint)  : Current derivative with respect to S of spatial
c                        rotation vectors at Gauss points for stiffness
c                        between
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: TSRF3F
c CALL     : SHP1D, TRCF3F, UPRM3F, UPOM3F, MSTS3F, BC3F, DGBG3F, PRSP3F
c            INT1D, MATQUA
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'
      include  'cdata.h'
      include  'counts.h'
      include  'ddata.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'erotas.h'
      include  'tdata.h'

      integer   ndf,ndm,nst,nh,isw, j1,k1,l1,m1,n1,k2,k3,m2
      real*8    dx,fac,xjac, temp, al0,al1, hn(*),h1(*)
      real*8    d(*),ul(ndf,nen,*),xl(ndm,nel),s(nst,nst),p(nst)
      real*8    omega(3,gint),thts(3,gint),thtpm(3,gint),xlbds(4,gint)
      real*8    sg(2,5),shp(2,3),exf(3),eforcl(3)
      real*8    mstrn(6),mstrs(6),sstrs(6),strs1(6)
      real*8    xnode(3,3),xnode0(3,3),xnode1(3,3),omg1(3),xlbd1(4)
      real*8    fi(3),fi0(3),fi1(3),fiprim(3),fi0prim(3),fi1prim(3)
      real*8    rotmat(3,3),rotmat0(3,3),rotmat1(3,3),rotmat2(3,3)
      real*8    aux1(3,3),aux2(3,3),aux3(3,3),aux4(3,3),xlbda0(3,3)
      real*8    aux5(3,3),aux6(3,3),aux7(3,3),aux8(3,3)
      real*8    btanc(6,18),btanc1(6,18),btanc2(6,18),bts1(6),bts2(6)
      real*8    shook(6,6),shook2(6,6),s1(12,12)
      real*8    bgdgbg(18,18),sfol(18,18)

      save

c     Recover current total translational displacements at t=n+alpha

      fac = 1.d0/theta(3) - 1.d0
      do j1 = 1,nel
        do k1 = 1,3
          xnode(k1,j1) = xl(k1,j1) + ul(k1,j1,1)
          if(theta(1).ne.0.d0) then
            xnode0(k1,j1) = xnode(k1,j1) - ul(k1,j1,2)
            xnode1(k1,j1) = xnode(k1,j1) + ul(k1,j1,2)*fac
          endif
        end do ! k1
      end do ! j1

c     Set assembly array

      m1 = 0
      n1 = 0
      do j1 = 1,nel
        do k1 = 1,6
          is(n1+k1) = m1 + k1
        end do ! k1
        m1 = m1 + ndf
        n1 = n1 + 6
      end do ! j1

c     Computation of Gauss points

      if(nint(d(182)).gt.0) then
        call int1dn(gint, sg)
      else
        call int1d(gint, sg)
      endif

c     Loop on each Gauss point

      do j1 = 1,gint
        do k1 = 1,12
          do k2 = 1,12
            s1(k2,k1) = 0.0d0
          end do ! k2
        end do ! k1
        call shp1d(sg(1,j1),xl,shp,ndm,nel,xjac)
        dx = sg(2,j1)*xjac

c       Residual forces and convected part of tangent stiffness matrix

c       Computation of current translational configuration and
c       its derivative with respect to S

        call trcf3f(xnode,nel,shp,fi,fiprim)
        if(theta(1).ne.0.d0) then
          call trcf3f(xnode0,nel,shp,fi0,fi0prim)
          call trcf3f(xnode1,nel,shp,fi1,fi1prim)
        endif

c       Computation of rotation matrix and curvature at t=n+alpha

        call uprm3f(xlbds(1,j1),thts(1,j1),xlbd1)
        call upom3f(omega(1,j1),thts(1,j1),thtpm(1,j1),omg1)

c       Computation of material strains & internal stresses at t=n+alpha
c       and spatial Hookean tensor

        call quamat(xlbds(1,j1),rotmat0)
        call quamat(xlbd1      ,rotmat1)

        al1 = theta(3)
        al0 = 1.d0 - al1
        do k1 =1,3
          do k2 =1,3
            if(theta(1).ne.0.d0) then
              rotmat(k1,k2)  = al0*rotmat0(k1,k2) + al1*rotmat1(k1,k2)
              rotmat2(k1,k2) = rotmat(k1,k2)
            else
              rotmat(k1,k2)  = rotmat1(k1,k2)
              rotmat2(k1,k2) = al0*rotmat0(k1,k2) + al1*rotmat1(k1,k2)
            endif
          end do ! k2
        end do ! k1

        if(theta(1).ne.0.d0) then
          call msts3f(d,fi0prim,fi1prim,rotmat0,rotmat1,omega(1,j1)
     &               ,omg1,mstrn,mstrs,hn,h1,nh,isw)
        else
          call msts3f(d,fiprim,fiprim,rotmat1,rotmat1,omg1
     &               ,omg1,mstrn,mstrs,hn,h1,nh,isw)
        endif

c       Save material stress for tplot use

        l1 = 6*(j1 -1)
        do k1 = 1,6
          tt(k1+l1) = mstrs(k1)
        end do !k1

c       Compute moduli for current state

        do k1 = 1,6
          do l1 = 1,6
            shook (l1,k1) = 0.0d0
            shook2(l1,k1) = 0.0d0
          end do ! l1
        end do ! k1

        do k1 = 1,3
          do l1 = 1,3
            aux2(k1,l1) = 0.0d0
            aux5(k1,l1) = 0.0d0
            aux4(k1,l1) = 0.0d0
            aux6(k1,l1) = 0.0d0
            do m1=1,3
              aux2(k1,l1) = aux2(k1,l1) + rotmat(k1,m1)*mhook(m1,l1)
              aux5(k1,l1) = aux5(k1,l1) + rotmat(k1,m1)*mhook(m1,l1+3)
              aux4(k1,l1) = aux4(k1,l1) + rotmat(k1,m1)*mhook(m1+3,l1+3)
              aux6(k1,l1) = aux6(k1,l1) + rotmat(k1,m1)*mhook(m1+3,l1)
            end do ! m1
          end do ! l1
        end do ! k1

        do k1 = 1,3
          do l1 = 1,3
            aux1(k1,l1) = 0.0d0
            aux3(k1,l1) = 0.0d0
            aux7(k1,l1) = 0.0d0
            aux8(k1,l1) = 0.0d0
            do m1=1,3
              shook(k1  ,l1) = shook(k1  ,l1)+aux2(k1,m1)*rotmat1(l1,m1)
              shook(k1+3,l1) = shook(k1+3,l1)+aux6(k1,m1)*rotmat1(l1,m1)
              aux1(k1,l1)  = aux1(k1,l1)  + aux4(k1,m1)*rotmat2(l1,m1)
              aux7(k1,l1)  = aux7(k1,l1)  + aux5(k1,m1)*rotmat2(l1,m1)
              aux3(k1,l1)  = aux3(k1,l1)
     &                     + aux4(k1,m1)*rotmat1(l1,m1)*0.5d0
              aux8(k1,l1)  = aux8(k1,l1)
     &                     + aux5(k1,m1)*rotmat1(l1,m1)*0.5d0
            end do ! m1
          end do ! l1
        end do ! k1

        call td3f(thts(1,j1),aux2)
        do k1 = 1,3
          shook2(k1  ,4) = aux8(k1,2)*thtpm(3,j1)
     &                   - aux8(k1,3)*thtpm(2,j1)
          shook2(k1  ,5) = aux8(k1,3)*thtpm(1,j1)
     &                   - aux8(k1,1)*thtpm(3,j1)
          shook2(k1  ,6) = aux8(k1,1)*thtpm(2,j1)
     &                   - aux8(k1,2)*thtpm(1,j1)

          shook2(k1+3,4) = aux3(k1,2)*thtpm(3,j1)
     &                   - aux3(k1,3)*thtpm(2,j1)
          shook2(k1+3,5) = aux3(k1,3)*thtpm(1,j1)
     &                   - aux3(k1,1)*thtpm(3,j1)
          shook2(k1+3,6) = aux3(k1,1)*thtpm(2,j1)
     &                   - aux3(k1,2)*thtpm(1,j1)

          shook (k1  ,4) = aux7(k1,1)*aux2(1,1)
     &                   + aux7(k1,2)*aux2(2,1)
     &                   + aux7(k1,3)*aux2(3,1)
          shook (k1  ,5) = aux7(k1,1)*aux2(1,2)
     &                   + aux7(k1,2)*aux2(2,2)
     &                   + aux7(k1,3)*aux2(3,2)
          shook (k1  ,6) = aux7(k1,1)*aux2(1,3)
     &                   + aux7(k1,2)*aux2(2,3)
     &                   + aux7(k1,3)*aux2(3,3)

          shook (k1+3,4) = aux1(k1,1)*aux2(1,1)
     &                   + aux1(k1,2)*aux2(2,1)
     &                   + aux1(k1,3)*aux2(3,1)
          shook (k1+3,5) = aux1(k1,1)*aux2(1,2)
     &                   + aux1(k1,2)*aux2(2,2)
     &                   + aux1(k1,3)*aux2(3,2)
          shook (k1+3,6) = aux1(k1,1)*aux2(1,3)
     &                   + aux1(k1,2)*aux2(2,3)
     &                   + aux1(k1,3)*aux2(3,3)
        end do ! k1

        call invert(aux2,3,3)

c       Spatial stresses and matrix Bc

        do k1 = 1,3
          sstrs(k1  ) = rotmat(k1,1)*mstrs(1)
     &                + rotmat(k1,2)*mstrs(2)
     &                + rotmat(k1,3)*mstrs(3)
          sstrs(k1+3) = rotmat(k1,1)*mstrs(4)
     &                + rotmat(k1,2)*mstrs(5)
     &                + rotmat(k1,3)*mstrs(6)
        end do ! k1
        call bc3f(shp,nel,fiprim,btanc)
        call bd3f(shp,nel,btanc2)
        if(theta(1).ne.0.d0)then
          do k1 = 1,3
            strs1(k1  ) = rotmat1(k1,1)*mstrs(1)
     &                  + rotmat1(k1,2)*mstrs(2)
     &                  + rotmat1(k1,3)*mstrs(3)
            strs1(k1+3) = rotmat1(k1,1)*mstrs(4)
     &                  + rotmat1(k1,2)*mstrs(5)
     &                  + rotmat1(k1,3)*mstrs(6)
          end do ! k1
          call bc3f(shp,nel,fi1prim,btanc1)
        endif

c       Computation of current external forces vector

        if(iforc.eq.1) then
          do k1 = 1,3
            exf(k1) = eforce(k1)
          end do ! k1
        else if(iforc.eq.2) then
          call quamat(rotqu0(1,j1), xlbda0)
          do k1 = 1,3
            exf(k1) = xlbda0(k1,1)*eforce(1)
     &              + xlbda0(k1,2)*eforce(2)
     &              + xlbda0(k1,3)*eforce(3)
          end do ! k1
        else if(iforc.eq.3) then
          call quamat(rotqu0(1,j1), xlbda0)
          do k1 = 1,3
            eforcl(k1) = xlbda0(1,k1)*eforce(1)
     &                 + xlbda0(2,k1)*eforce(2)
     &                 + xlbda0(3,k1)*eforce(3)
          end do ! k1
          do k1 = 1,3
            exf(k1)    = rotmat(k1,1)*eforcl(1)
     &                 + rotmat(k1,2)*eforcl(2)
     &                 + rotmat(k1,3)*eforcl(3)
          end do ! k1
        else if(iforc.eq.4) then
          do k1 = 1,3
            exf(k1)    = rotmat(k1,1)*eforce(1)
     &                 + rotmat(k1,2)*eforce(2)
     &                 + rotmat(k1,3)*eforce(3)
          end do ! k1
        end if

c       Distributed applied loading

        k1 = 0
        do l1 = 1,nel
          temp = shp(2,l1)*dx
          p(k1+1) = p(k1+1) + exf(1)*temp
          p(k1+2) = p(k1+2) + exf(2)*temp
          p(k1+3) = p(k1+3) + exf(3)*temp
          k1      = k1 + ndf
        end do ! l1

c       Gravity loading

        call fbody3d(d,xl, p, ndm,ndf, isw)

c       Subtract stress divergence term

        do l1 = 1,6
          temp = sstrs(l1)*dx
          do k1 = 1,6*nel
            p(k1) = p(k1) - btanc(l1,k1)*temp
          end do ! k1
        end do ! l1

        do k1 = 1,6*nel
          do l1 = 1,6
            bts1(l1) = 0.0d0
            bts2(l1) = 0.0d0
            do m1 = 1,6
              bts1(l1) = bts1(l1) + btanc(m1,k1)*shook (m1,l1)*dx
              bts2(l1) = bts2(l1) + btanc(m1,k1)*shook2(m1,l1)*dx
            end do ! m1
          end do ! l1

          if(theta(1).ne.0.d0)then
            do n1 = 1,6*nel
              do l1 = 1,6
                s1(k1,n1) = s1(k1,n1) + bts1(l1)*btanc1(l1,n1)
     &                                + bts2(l1)*btanc2(l1,n1)
              end do ! n1
            end do ! l1
          else
            do n1 = 1,6*nel
              do l1 = 1,6
                s1(k1,n1) = s1(k1,n1) + bts1(l1)*btanc (l1,n1)
     &                                + bts2(l1)*btanc2(l1,n1)
              end do ! n1
            end do ! l1
          endif
        end do ! k1

c       Geometric part of tangent stiffness matrix

c         Computation of matrix Bgt*Dg*Bg

          if(theta(1).ne.0.d0)then
            call dgbg3f(shp,nel,fiprim,sstrs,strs1,bgdgbg)
          else
            call dgbg3f(shp,nel,fiprim,sstrs,sstrs,bgdgbg)
          endif

c         Geometric part of tangent stiffness matrix

          do k1 = 1,6*nel
            do l1 = 1,6*nel
              s1(k1,l1) = s1(k1,l1) + bgdgbg(k1,l1)*dx
            end do ! l1
          end do ! k1
c       end if

c       Follower forces components of tangent stiffness matrix

        if(iforc.eq.3.or.iforc.eq.4) then
          if(iforc.eq.3) then
            call quamat(rotqu0(1,j1), xlbda0)
            do k1 = 1,3
              eforcl(k1) = xlbda0(1,k1)*eforce(1)
     &                   + xlbda0(2,k1)*eforce(2)
     &                   + xlbda0(3,k1)*eforce(3)
            end do ! k1
            do k1 = 1,3
              exf(k1)    = rotmat1(k1,1)*eforcl(1)
     &                   + rotmat1(k1,2)*eforcl(2)
     &                   + rotmat1(k1,3)*eforcl(3)
            end do ! k1
          else if(iforc.eq.4) then
            do k1 = 1,3
              exf(k1)    = rotmat1(k1,1)*eforce(1)
     &                   + rotmat1(k1,2)*eforce(2)
     &                   + rotmat1(k1,3)*eforce(3)
            end do ! k1
            call fol13f(shp,nel,exf,sfol)
          end if
          do k1 = 1,6*nel
            do l1 = 1,6*nel
              s1(k1,l1) = s1(k1,l1) + sfol(k1,l1)*dx
            end do ! l1
          end do ! k1
        end if

        m1 = 0
        do m2 = 1,nel
          call td3f(rots(1,m2,2),aux1)
          do k2 = 1,3
            do k1 = 1,3
              aux3(k1,k2) = aux2(k1,1)*aux1(1,k2)
     &                    + aux2(k1,2)*aux1(2,k2)
     &                    + aux2(k1,3)*aux1(3,k2)
            end do ! k1
          end do ! k2

          do k1 = 1,6*nel
            do k2 = 1,3
              s(is(k1),is(k2+m1)) = s(is(k1),is(k2+m1)) + s1(k1,k2+m1)
              do k3 = 1,3
                s(is(k1),is(k2+m1+3)) = s(is(k1),is(k2+m1+3))
     &                           + s1(k1,k3+m1+3)*aux3(k3,k2)
              end do ! k3
            end do ! k2
          end do ! k1
          m1 = m1 + 6
        end do ! m2
      end do ! j1

c     Compute rotational body force residual and tangent

      if(d(4).gt.0.0d0 .and. d(65).gt.0.0d0) then
        call fbody3w(d(4)*d(32),d(65),xl,ul, p,s, isw.eq.3)
      endif

c     Multiply tangent by time integration factor

      if(theta(1).ne.0.d0) then
        do k1 = 1,ndf*nel
          do k2 = 1,ndf*nel
            s(k1,k2) = s(k1,k2)*ctan(1)
          end do ! k2
        end do ! k1
      endif

      end

      subroutine dyna3f(ul,xl,s,p,ndf,ndm,nst)

c     Computes dynamic part of tangent stiffness and residual forces
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: TSRF3F
c CALL     : SHP1D, UPRM3F, TD3F, PRSP3F
c            INT1D, INT1DL, MATQUA
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'
      include  'cdata.h'
      include  'ddata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'erotas.h'
      include  'tdata.h'

      integer   ndf,ndm,nst
      integer   j1,k1,l1,m1,k2,l2,m2
      real*8    ul(ndf,nen,*),xl(ndm,nel),s(nst,nst),p(ndf,*)

      real*8    sg(2,5),shp(2,3),xjac,dx,t3(3),aa,dx1,dx2,vk,vn
      real*8    sinn1(3,3,3),sinn(3,3,3),wsnn(3,3),wsnn1(3,3),acc(3)
      real*8    rotnn(3,3,3),rotnn1(3,3,3),rotm(3,3,3)
      real*8    angnn1(3,3),angnn(3,3)
      real*8    aux1(3,3),aux2(3,3),aux3(3,3),aux4(3,3)

      save

c     Transfer of spatial velocities at nodes at t=n+1 & computation of
c     velocities at t=n

      do j1 = 1,nel

        do k1 = 1,3
          wsnn1(k1,j1) = rvel(k1,j1,2)   ! At t=n+1
          wsnn (k1,j1) = rvel(k1,j1,1)   ! At t=n
        end do ! k1

c       Computation of rotation matrices of nodes at t=n and t=n+1

        call quamat(xln(1,j1,1),rotnn(1,1,j1))
        call quamat(xln(1,j1,3),rotnn1(1,1,j1))

c       Computation of inertia dyadic at nodes at t=n and t=n+1

        do l1 = 1,3
          do k1 = 1,3
            aux1(k1,l1) = rotnn1(k1,1,j1)*miner(1,l1)
     &                  + rotnn1(k1,2,j1)*miner(2,l1)
     &                  + rotnn1(k1,3,j1)*miner(3,l1)

            aux2(k1,l1) = rotnn (k1,1,j1)*miner(1,l1)
     &                  + rotnn (k1,2,j1)*miner(2,l1)
     &                  + rotnn (k1,3,j1)*miner(3,l1)
          end do ! k1
        end do ! l1

        do l1 = 1,3
          do k1 = 1,3
            sinn1(k1,l1,j1) = aux1(k1,1)*rotnn1(l1,1,j1)
     &                      + aux1(k1,2)*rotnn1(l1,2,j1)
     &                      + aux1(k1,3)*rotnn1(l1,3,j1)

            sinn (k1,l1,j1) = aux2(k1,1)*rotnn (l1,1,j1)
     &                      + aux2(k1,2)*rotnn (l1,2,j1)
     &                      + aux2(k1,3)*rotnn (l1,3,j1)

            rotm (k1,l1,j1) = rotnn1(1,k1,j1)*rotnn(1,l1,j1)
     &                      + rotnn1(2,k1,j1)*rotnn(2,l1,j1)
     &                      + rotnn1(3,k1,j1)*rotnn(3,l1,j1)
          end do ! k1
        end do ! l1

c       Computation of angular momentum at t=n+1

        do k1 = 1,3
          angnn1(k1,j1) = sinn1(k1,1,j1)*wsnn1(1,j1)
     &                  + sinn1(k1,2,j1)*wsnn1(2,j1)
     &                  + sinn1(k1,3,j1)*wsnn1(3,j1)
          angnn (k1,j1) = sinn (k1,1,j1)*wsnn (1,j1)
     &                  + sinn (k1,2,j1)*wsnn (2,j1)
     &                  + sinn (k1,3,j1)*wsnn (3,j1)
        end do ! k1

      end do ! j1

c     Gauss points

      if(lobatt.eq.0) then
        call int1d (mint,sg)
      else
        call int1dl(mint,sg)
      end if

c     Loop over Gauss points for mass

      aa = theta(2)/(theta(1)*dt)
      do j1 = 1,mint

c       Shape functions

        call shp1d(sg(1,j1),xl,shp,ndm,nel,xjac)
        xjac = sg(2,j1)*xjac
        dx   = xjac/(theta(2)*dt)
        dx1  = arho*xjac
        dx2  = ctan(3)*dx1

c       Residual: Compute mass * acceleration at integration point

        do l1 = 1,3
          acc(l1) = 0.0d0
          do m1 = 1,nel
            acc(l1) = acc(l1) + shp(2,m1)*ul(l1,m1,5)
          end do ! m1
          acc(l1) = acc(l1)*dx1
        end do ! l1

        do k1 = 1,nel
          vk = dx*shp(2,k1)
          do l1 = 1,3
            p(l1,k1)   = p(l1,k1)   - shp(2,k1)*acc(l1)
            p(l1+3,k1) = p(l1+3,k1) - vk*(angnn1(l1,k1)-angnn(l1,k1))
          end do ! l1
        end do ! k1

c       Computation of dynamic part of tangent stiffness matrix

c       Translational part

        k2 = 0
        do k1 = 1,nel
          vk = dx2*shp(2,k1)
          l2 = 0
          do l1 = 1,nel
            vn = vk*shp(2,l1)
            do m1 = 1,3
              s(k2+m1,l2+m1) = s(k2+m1,l2+m1) + vn
            end do ! m1
            l2 = l2 + ndf
          end do ! l1
          k2 = k2 + ndf
        end do ! k1

c       Rotational part

        l2 = 3
        do l1 = 1,nel
          call hat3f(angnn1(1,l1),aux1)
          call hat3f(wsnn1(1,l1),aux4)

          do k1 = 1,3
            t3(k1) = rotm (k1,1,l1)*wsnn (1,l1)
     &             + rotm (k1,2,l1)*wsnn (2,l1)
     &             + rotm (k1,3,l1)*wsnn (3,l1)
            do k2 = 1,3
              aux2(k1,k2) = sinn1(k1,1,l1)*aux4(1,k2)
     &                    + sinn1(k1,2,l1)*aux4(2,k2)
     &                    + sinn1(k1,3,l1)*aux4(3,k2)
            end do ! k2
          end do ! k1

          call hat3f(t3,aux4)
          call td3f(rots(1,l1,2),aux3)

          do k1 = 1,3
            do k2 = 1,3
              aux4(k1,k2) = aux4(k1,k2) + aa*aux3(k1,k2)
            end do ! k2
          end do ! k1

          do k1 = 1,3
            do k2 = 1,3
              aux3(k1,k2) = sinn1(k1,1,l1)*aux4(1,k2)
     &                    + sinn1(k1,2,l1)*aux4(2,k2)
     &                    + sinn1(k1,3,l1)*aux4(3,k2)
            end do ! k2
          end do ! k1

          vk = dx*shp(2,l1)
          do m1 = 1,3
            do m2 = 1,3
              s(m1+l2,m2+l2) = s(m1+l2,m2+l2) + vk*(-aux1(m1,m2)
     &                          + aux2(m1,m2) + aux3(m1,m2))
            end do ! m2
          end do ! m1
          l2 = l2 + ndf
        end do ! l1
      end do ! j1

      end

      subroutine trcf3f(xnode,nel,shp,fi,fiprim)

c     Computes current position of one point of line of centroids and
c     its derivative with respect to S at that point
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i     xnode(3,nel)   : Translational coordinate of nodes
c   o   fi(3)          : Current position of line of centroids
c   o   fiprim(3)      : Derivative of fi with respect to S at point
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: STAT3F, OUTP3F
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nel,j1,k1
      real*8    xnode(3,nel),shp(2,nel),fi(3),fiprim(3)

c     Computation of fi and fiprim

      do k1 = 1,3
        fi(k1)     = xnode(k1,1)*shp(2,1) + xnode(k1,2)*shp(2,2)
        fiprim(k1) = xnode(k1,1)*shp(1,1) + xnode(k1,2)*shp(1,2)
        do j1 = 3,nel
          fi(k1)     = fi(k1)     + shp(2,j1)*xnode(k1,j1)
          fiprim(k1) = fiprim(k1) + shp(1,j1)*xnode(k1,j1)
        end do ! j1
      end do ! k1

      end

      subroutine msts3f(d,fi0prim,fi1prim,rotmat0,rotmat1,omg0
     &                  ,omg1,strain,stress,hn,h1,nh,isw)

c      Computes material strains and stresses at one point
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i     fi0prim(3)      : Derivative of fi with respect to S at t_n
c i     fi1prim(3)      : Derivative of fi with respect to S at t_n+1
c i     rotmat0(9)      : Rotation matrix at that point at t_n
c i     rotmat1(9)      : Rotation matrix at that point at t_n+1
c i     omg0(3)         : Spatial curvature at that point at t_n
c i     omg1(3)         : Spatial curvature at that point at t_n+1
c   o   strain(6)      : Material strains
c   o   stress(6)      : Material stresses
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: STAT3F, OUTP3F
c CALL     : BM3RES
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'
      include  'ddata.h'

      integer   i, j,isw, nh
      real*8    al0,al1, d(*),hn(*),h1(*),omg0(3),omg1(3)
      real*8    fi0prim(3),fi1prim(3),rotmat0(3,3),rotmat1(3,3)
      real*8    strain0(6),strain1(6),strain(6),stress(6)

      save

c     Computation of strains

      do i = 1,3
        strain0(i  ) = 0.0d0
        strain0(i+3) = 0.0d0
        strain1(i  ) = 0.0d0
        strain1(i+3) = 0.0d0
        do j = 1,3
          strain0(i  ) = strain0(i  ) + rotmat0(j,i)*fi0prim(j)
          strain0(i+3) = strain0(i+3) + rotmat0(j,i)*omg0(j)
          strain1(i  ) = strain1(i  ) + rotmat1(j,i)*fi1prim(j)
          strain1(i+3) = strain1(i+3) + rotmat1(j,i)*omg1(j)
        end do ! j
      end do ! j

      al1 = theta(3)
      al0 = 1.d0 - al1
      do i = 1,6
        strain(i) = al0 * strain0(i) + al1 * strain1(i)
      end do ! i
      strain(3) = strain(3) - 1.d0

c     Computation of stresses

      do i = 1,6
        stress(i) = 0.0d0
        do j = 1,6
          stress(i) = stress(i) + mhook(i,j)*strain(j)
        end do ! j
      end do ! i

      if(d(101).gt.0.0d0) then
        call bm3res(d,hn,h1,nh,strain, stress,mhook, isw)
      endif

      end

      subroutine bc3f(shp,nel,fiprim,btanc)

c     Computes matrix Bt corresponding to convected part of tangent
c     stiffness matrix at a certain Gauss point
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i     fiprim(3)      : Derivative of fi with respect to S at point
c   o   btanc(18,6)    : Matrix Bt
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: STAT3F
c CALL     : VECP
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nel,j1,k1,j2,k2
      real*8    shp(2,3),fiprim(3),btanc(6,18),t1(3)

      save

c     Computation of btanc

      do j1 = 1,18
        do k1 = 1,6
          btanc(k1,j1) = 0.0d0
        end do ! k1
      end do ! j1

      j2 = 0
      do j1 = 1,nel
        do k1 = 1,3
          k2 = k1 + j2
          btanc(k1  ,k2  ) = shp(1,j1)
          btanc(k1+3,k2+3) = shp(1,j1)
          t1( 1) = 0.0d0
          t1( 2) = 0.0d0
          t1( 3) = 0.0d0
          t1(k1) = shp(2,j1)
          call vecp(fiprim,t1,btanc(1,k2+3))
        end do ! k1
        j2 = j2 + 6
      end do ! j1

      end

      subroutine bd3f(shp,nel,btanc)

c     Computes matrix Bt corresponding to convected part of tangent
c     stiffness matrix at a certain Gauss point
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c   o   btanc(18,6)    : Matrix Bd
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: STAT3F
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nel,j1,k1,j2
      real*8    shp(2,3),btanc(6,18)

c     Computation of btanc

      do j1 = 1,6*nel
        do k1 = 1,6
          btanc(k1,j1) = 0.0d0
        end do ! k1
      end do ! j1
      j2 = 0
      do j1 = 1,nel
        do k1 = 1,3
          btanc(k1+3,k1+j2+3) = shp(2,j1)
        end do ! k1
        j2 = j2 + 6
      end do ! j1

      end

      subroutine dgbg3f(shp,nel,fiprim,strs,strs1,bgdgbg)

c     Computes product BgT*Dg*Bg corresponding to geometric part of
c     tangent stiffness matrix
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i     fiprim(3)      : Derivative of fi with respect to S at point
c i     strs(6)        : Spatial internal stresses
c i     strs1(6)       : modified Spatial internal stresses
c   o   bgdgbg(9,9)    : Matrix BgT*Dg*Bg
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: STAT3F
c CALL     : BG3F, MULTFD, VECP
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nel,j1,k1,j2,k2,l1,l2,m1
      real*8    shp(2,3),fiprim(3),strs(6),strs1(6),bgdgbg(18,18)
      real*8    btang(162),btangt(162),dgbg(162),aux(3,3),scalar

      save

c     Computation of matrix Bg and its transpose

      call bg3f(shp,nel,btang,btangt)

c     Computation of formal product Dg*Bg

c     Scalar and tensor product of n*fip

      do j1 = 1,162
        dgbg(j1) = 0.0d0
      end do ! j1
      scalar = strs1(1)*fiprim(1)+strs1(2)*fiprim(2)+strs1(3)*fiprim(3)
      do j1 = 1,3
        do k1 = 1,3
          aux(j1,k1) = strs1(j1)*fiprim(k1)
        end do ! k1
        aux(j1,j1) = aux(j1,j1) - scalar
      end do ! j1

c     Matrix Dg*Bg

      j2 = 1
      do j1 = 1,nel
        do k1 = 27,45,9
          k2 = k1 + j2
          call vecp(btang(k2+6),strs1   ,dgbg(k2  ))
          call vecp(btang(k2+6),strs1(4),dgbg(k2+3))
          call vecp(strs   ,btang(k2-27),dgbg(k2-21))
          do l1 = 1,3
            l2 = l1 + k2 + 5
            do m1 = 1,3
              dgbg(l2) = dgbg(l2) + aux(l1,m1)*btang(m1+k2+5)
            end do ! m1
          end do ! l1
        end do ! k1
        j2 = j2 + 54
      end do ! j1

c     Matrix BgT*Dg*Bg

      call multfd(btangt,dgbg,bgdgbg,18,9,18)

      end

      subroutine bg3f(shp,nel,btang,btangt)

c     Computes matrix Bg corresponding to geometric part of
c     tangent stiffness matrix and its transpose
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c   o   btang(18,9)    : Matrix Bg
c   o   btangt(9,18)   : Matrix BgT
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: DGBG3F
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nel,j1,j2,k1
      real*8    shp(2,3),btang(9,18),btangt(18,9)

c     Computation of Bg

      do j1 = 1,6*nel
        do k1 = 1,9
          btang(k1,j1) = 0.0d0
        end do ! k1
      end do ! j1
      j2 = 0
      do j1 = 1,nel
        do k1 = 1,3
          btang(k1  ,j2+k1  ) = shp(1,j1)
          btang(k1+3,j2+k1+3) = shp(1,j1)
          btang(k1+6,j2+k1+3) = shp(2,j1)
        end do ! k1
        j2 = j2 + 6
      end do ! j1

c     Computation of BgT

      do j1 = 1,9
        do k1 = 1,18
          btangt(k1,j1) = btang(j1,k1)
        end do ! k1
      end do ! j1

      end

      subroutine mass3f(xl,s,p,ndf,ndm,nst)

c     Computes mass matrix
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: TSRF3F
c CALL     : SHP1D, UPRM3F, TD3F, PRSP3F
c            INT1D, INT1DL, MATQUA
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'
      include  'eldata.h'
      include  'erotas.h'

      integer   ndf,ndm,nst,j1,k1,l1,m1,k2,l2,m2
      real*8    xjac,dx1,vm,vn,xl(ndm,nel),s(nst,nst),p(ndf,*)
      real*8    sg(2,5),shp(2,3),aux(3,3)
      real*8    rotnn(3,3,3),sinn(3,3,3),tinn(3,3)

      save

      do j1 = 1,nel

c       Computation of rotation matrices of nodes at t=n+1

        call quamat(xln(1,j1,3),rotnn(1,1,j1))

c       Computation of inertia dyadic at nodes at t=n+1

        do l1 = 1,3
          do k1 = 1,3
            aux(k1,l1) = rotnn(k1,1,j1)*miner(1,l1)
     &                 + rotnn(k1,2,j1)*miner(2,l1)
     &                 + rotnn(k1,3,j1)*miner(3,l1)
          end do ! k1
        end do ! l1

        do l1 = 1,3
          do k1 = 1,3
            tinn(k1,l1) = aux(k1,1)*rotnn(l1,1,j1)
     &                  + aux(k1,2)*rotnn(l1,2,j1)
     &                  + aux(k1,3)*rotnn(l1,3,j1)
          end do ! k1
        end do ! l1

        call td3f(rots(1,j1,2),aux)

        do l1 = 1,3
          do k1 = 1,3
            sinn(k1,l1,j1) = tinn(k1,1)*aux(1,l1)
     &                     + tinn(k1,2)*aux(2,l1)
     &                     + tinn(k1,3)*aux(3,l1)
          end do ! k1
        end do ! l1

      end do ! j1

c     Gauss points

      if(lobatt.eq.0) then
        call int1d (mint,sg)
      else
        call int1dl(mint,sg)
      end if

c     Loop over Gauss points for mass

      do j1 = 1,mint

c       Shape functions

        call shp1d(sg(1,j1),xl,shp,ndm,nel,xjac)
        xjac = sg(2,j1)*xjac
        dx1  = arho*xjac

c       Computation of mass matrix

c       Translational part

        k2 = 0
        do k1 = 1,nel
          vm = dx1*shp(2,k1)

c         Lumped mass part

          do m1 = 1,3
            p(m1,k1) = p(m1,k1) + vm
          end do ! m1

c         Consistent mass part

          l2 = 0
          do l1 = 1,nel

            vn = vm*shp(2,l1)

            do m1 = 1,3
              s(k2+m1,l2+m1) = s(k2+m1,l2+m1) + vn
            end do ! m1
            l2 = l2 + ndf

          end do ! l1
          k2 = k2 + ndf

        end do ! k1

c       Rotational part

        l2 = 3
        do l1 = 1,nel

          vm = xjac*shp(2,l1)

          do m1 = 1,3
            do m2 = 1,3
              s(l2+m1,l2+m2) = s(l2+m1,l2+m2) + vm*sinn(m1,m2,l1)
            end do ! m2
          end do ! m1
          l2 = l2 + ndf

        end do ! l1
      end do ! j1

      end

      subroutine fol13f(shp,nel,exf,sfol)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Computes part of tangent stiffness matrix corresponding
c     to follower forces along directors

c - ARGUMENTS
c i     exf(3)         : Current external forces vector in global coord
c   o   sfol(nst,nst)  : Tangent stiffness matrix (part coresponding to
c                        follower forces
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: STAT3F
c CALL     :
c            PZERO, VECP
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nel,j1,k1,j2,k2,l1,l2,m1
      real*8    shp(2,3),exf(3),sfol(18,18),t1(3),t2(3)

      save

c     Computation of follower forces part of tangent stiffness matrix

      do j1 = 1,6*nel
        do k1 = 1,6*nel
          sfol(k1,j1) = 0.0d0
        end do ! k1
      end do ! j1
      j2 = 3
      do j1 = 1,nel
        do k1 =  1,3
          k2 = j2 + k1
          t1( 1) = 0.0d0
          t1( 2) = 0.0d0
          t1( 3) = 0.0d0
          t1(k1) = shp(2,j1)
          call vecp(exf,t1,t2)
          l2 = 0
          do l1 = 1,nel
            do m1 = 1,3
              sfol(m1+l2,k2) = shp(2,l1)*t2(m1)
            end do ! m1
            l2 = l2 + 6
          end do ! l1
        end do ! k1
        j2 = j2 + 6
      end do ! j1

      end

      subroutine fol23f(shp,nel,eforce,centr,rmod3,rmod2,rotmat,sfol)

c     Computes part of tangent stiffness matrix corresponding to
c     follower forces along line of centroids (external preassure)

c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i     eforce(3)      : External forces vector in current coordinates
c i     centr(3,3)     : Current orthonormal frame along which load is
c                        defined
c i     rmod3,rmod2    : Moduli of nonunitary vectors 3 and 2 of centr
c i     rotmat(3,3)    : The director frame
c   o   sfol(nst,nst)  : Tangent stiffness matrix (part coresponding to
c                        follower forces
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: STAT3F
c CALL     :
c            PZERO, VECP, MULTFD
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nel,j1,k1,j2,k2,k3,l1,l2,m1
      real*8    shp(2,3),eforce(3),centr(3,3),rmod3,rmod2,rotmat(3,3)
      real*8    sfol(18,18),t1(3),t2(3),t3(3),t4(3),da(3,3),scal1,scal2

      save

c     Computation of follower forces part of tangent stiffness matrix

      do j1 = 1,6*nel
        do k1 = 1,6*nel
          sfol(k1,j1) = 0.0d0
        end do ! k1
      end do ! j1
      j2 = 0
      do j1 = 1,nel
        do k1 = 1,6
          k2 = j2 + k1
          k3 = k1 - 3*(k1/4)
          do l1 = 1,3
            t1(l1) = 0.0d0
          end do ! l1
          t1(k3) = shp(2,j1)

c         Computations of variations of three vectors composing centr

c         The tangent vector to line of centroids (a3)

          if(k1.le.3) then
            do l1 = 1,3
              da(l1,3) = (t1(l1) - t1(k3)*centr(k3,3)*centr(l1,3))/rmod3
            end do ! l1
          else
            do l1 = 1,3
              da(l1,3) = 0.0d0
            end do ! l1
          end if

c         The second vector orthogonal to a3 (a2)

          scal1 = rotmat(1,2)*centr(1,3)
     &          + rotmat(2,2)*centr(2,3)
     &          + rotmat(3,2)*centr(3,3)
          scal2 = rotmat(1,2)*da(1,3)
     &          + rotmat(2,2)*da(2,3)
     &          + rotmat(3,2)*da(3,3)
          do l1 = 1,3
            t2(l1) = rmod2*centr(k1,2)
          end do ! l1
          call vecp(t1,t2,t4)
          do l1 = 1,3
            if(k1.le.3) then
              t3(l1) = -scal1*da(l1,3) - scal2*centr(l1,3)
            else
              t3(l1) = t4(l1) - t4(k3)*centr(k3,3)*centr(l1,3)
            end if
          end do ! l1
          scal2 = t2(1)*t3(1) + t2(2)*t3(2) + t2(3)*t3(3)
          do l1 = 1,3
            da(l1,2) = t3(l1)/rmod2 - scal2*t2(l1)/rmod2**3
          end do ! l1

c         The first vector orthogonal to a3 (a1)

          call vecp(da(1,2),centr(1,3),t2)
          call vecp(centr(1,2),da(1,3),t3)
          do l1 = 1,3
            da(l1,1) = t2(l1) + t3(l1)
          end do ! l1

c         Computation of total variation

          do m1 = 1,3
            t2(m1) = da(m1,1)*eforce(1)
     &             + da(m1,2)*eforce(2)
     &             + da(m1,3)*eforce(3)
          end do ! m1
          l2 = 0
          do l1 = 1,nel
            do m1 = 1,3
              sfol(m1+l2,k2) = shp(2,l1)*t2(m1)
            end do ! m1
            l2 = l2 + 6
          end do ! l1
        end do ! k1
        j2 = j2 + 6
      end do ! j1

      end

      subroutine indb3f(xlbds,xlbdd,omega,thts,thtpm,thtd)

c     Initializes database at start of whole process

c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c   o   xlbds(4,gint)  : Quaternions associated with rotation matrices
c                        at Gauss points for stiffness at t=n
c   o   xlbdd(4,mint)  : Quaternions associated with rotation matrices
c                        at Gauss points for mass at t=n
c   o   omega(3,gint)  : Curvatures at Gauss points of stiffness at t=n
c   o   thts(3,gint)   : Current spatial rotation vectors at Gauss
c                        points for stiffness between t=n and t=n+alpha
c   o   thtpm(3,gint)  : Current derivative with respect to S: spatial
c                        rotation vectors at Gauss points for stiffness
c                        between t=n and t=n+alpha
c   o   thtd(3,mint)   : Current spatial rotation vectors at Gauss
c                        points for mass between t=n and t=n+alpha
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: FRM3FD, TSRF3F
c CALL     : QUATEX, PZERO, PMOVE, PRCM3F
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'
      include  'eldata.h'
      include  'erotas.h'

      integer   j1,k1
      real*8    xlbds(4,gint),xlbdd(4,mint),omega(3,gint),thts(3,gint)
      real*8    thtpm(3,gint),thtd(3,mint)

c     Rotation matrix is initialized for all points of element

      do j1 = 1,gint
        do k1 = 1,4
          xlbds(k1,j1) = rotqu0(k1,j1)
        end do ! j1
        do k1 = 1,3
          omega(k1,j1) = 0.0d0
          thts(k1,j1)  = 0.0d0
          thtpm(k1,j1) = 0.0d0
        end do ! j1
      end do ! k1

      do j1 = 1,mint
        do k1 = 1,4
          xlbdd(k1,j1) = rotqu0(k1,j1+gint)
        end do ! j1
        do k1 = 1,3
          thtd(k1,j1) = 0.0d0
        end do ! j1
      end do ! k1

      end

      subroutine updb3f(d,xl,ndm,thts,thtpm)

c     Update database along iteration process
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i o   thts(3,gint)   : Current spatial rotation vectors at Gauss
c                        points for stiffness between t=n and t=n+alpha
c i o   thtpm(3,gint)  : Current derivative with respect to S of spatial
c                        rotation vectors at Gauss points for stiffness
c                        between t=n and t=n+alpha
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: FRM3FD, TSRF3F
c CALL     : SHP1D, UPRD3F, UPRV3F, PRCM3F
c            INT1D, INT1DL, PZERO, PMOVE,
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'
      include  'cdata.h'
      include  'ddata.h'
      include  'eldata.h'
      include  'erotas.h'
      include  'tdata.h'

      integer   ndm, j1,k1,k2
      real*8    d(*),xl(ndm,nel),thts(3,gint),thtpm(3,gint)
      real*8    sg(2,5),shp(2,3),xjac

      save

c     Update of rotation vectors at Gauss points for stiffness

c     Gauss points and pointers

      if(nint(d(182)).gt.0) then
        call int1dn(gint, sg)
      else
        call int1d(gint, sg)
      endif

c     Computation of incremental rotation vector and its derivative
c     with respect to S at t=n+1

      do j1 = 1,gint
        call shp1d(sg(1,j1),xl,shp,ndm,nel,xjac)
        do k2 = 1,3
          thts(k2,j1)  = 0.0d0
          thtpm(k2,j1) = 0.0d0
          do k1 = 1,nel
            thts (k2,j1) = thts (k2,j1) + rots(k2,k1,2)*shp(2,k1)
            thtpm(k2,j1) = thtpm(k2,j1) + rots(k2,k1,2)*shp(1,k1)
          end do ! k1
        end do ! k2
      end do ! j1

      end

      subroutine uprm3f(xlbda,tht,xlbdp)

c     Update rotation matrix given a certain incremental rotation vector
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i       xlbda(4) - Quaternions associated to initial rotation matrix
c i       tht(3)   - Spatial incremental rotation vector
c   o     xlbdp(4) - Quaternions associated to final rotation matrix
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: TSDB3F, STAT3F, DYNA3F
c CALL     : QRCAY, QUAMUL
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    xlbda(4),tht(3),xlbdp(4), qtht(4)

      save

c     Computation of quaternions associated to tht

      call rqcay(tht,qtht)

c     Computation of quaternions of new rotation matrix

      call quamul(qtht,xlbda,xlbdp)

      end

      subroutine upom3f(omega,rot,rotp,omegp)

c     Update spatial curvature given ertain incremental rotation vector
c     and its derivative with respect to S
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i       omega(3) - Spatial initial curvature
c i       tht(3)   - Spatial incremental rotation vector
c i       thtprm(3)- Derivative with respect to S of spatial incremental
c                    rotation vector
c   o     omegp(3) - Spatial final curvature
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: TSDB3F, STAT3F, DYNA3F
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    det,q31, q32, q33
      real*8    omega(3),omegp(3),rot(3),rotp(3),q(3,3),rr(3)

      rr(1) = rotp(1)+omega(1)+.5d0*(rot(2)*omega(3)-rot(3)*omega(2))
      rr(2) = rotp(2)+omega(2)+.5d0*(rot(3)*omega(1)-rot(1)*omega(3))
      rr(3) = rotp(3)+omega(3)+.5d0*(rot(1)*omega(2)-rot(2)*omega(1))

      q(1,1) =  1.0d0
      q(1,2) =  0.5d0*rot(3)
      q(1,3) = -0.5d0*rot(2)

      q(2,1) = -0.5d0*rot(3)
      q(2,2) =  1.0d0
      q(2,3) =  0.5d0*rot(1)

      q(3,1) =  0.5d0*rot(2)
      q(3,2) = -0.5d0*rot(1)
      q(3,3) =  1.0d0

      q31    = q(1,2)*q(2,3) - q(2,2)*q(1,3)
      q32    = q(1,3)*q(2,1) - q(2,3)*q(1,1)
      q33    = q(1,1)*q(2,2) - q(2,1)*q(1,2)

      det    = 1.d0/(q31*q(3,1) + q32*q(3,2) + q33*q(3,3))

      omegp(1) = (rr(1)*(q(2,2)*q(3,3) - q(3,2)*q(2,3))
     &         +  rr(2)*(q(3,2)*q(1,3) - q(1,2)*q(3,3))
     &         +  rr(3)* q31 )*det

      omegp(2) = (rr(1)*(q(2,3)*q(3,1) - q(3,3)*q(2,1))
     &         +  rr(2)*(q(3,3)*q(1,1) - q(1,3)*q(3,1))
     &         +  rr(3)* q32 )*det

      omegp(3) = (rr(1)*(q(2,1)*q(3,2) - q(3,1)*q(2,2))
     &         +  rr(2)*(q(3,1)*q(1,2) - q(1,1)*q(3,2))
     &         +  rr(3)* q33 )*det

      end

      subroutine tsdb3f(xlbds,xlbdd,omega,thts,thtpm,thtd)

c     Initializes database at start of a new time step
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i o   xlbda(4,gint)  : Quaternions associated with rotation matrices
c                        at Gauss points for stiffness at t=n
c i o   xlbdd(4,mint)  : Quaternions associated with rotation matrices
c                        at Gauss points for mass at t=n
c i o   omega(3,gint)  : Curvatures at Gauss points of stiffness at t=n
c i o   thts(3,gint)   : Current spatial rotation vectors at Gauss
c                        points for stiffness between t=n and t=n+alpha
c i o   thtpm(3,gint)  : Current derivative of spatial rotation vectors
c                        at Gauss points for stiffness between t=n and
c                        t=n+alpha
c i o   thtd(3,mint)   : Current spatial rotation vectors at Gauss
c                        points for mass between t=n and t=n+alpha
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: FRM3FD
c CALL     : UPRM3F, UPOM3F
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'
      include  'comblk.h'
      include  'cdata.h'
      include  'ddata.h'
      include  'eldata.h'
      include  'hdata.h'
      include  'tdata.h'

      integer   ilbd,iom,j1,k1
      real*8    t1(3),t2(3),newdb(100)
      real*8    xlbds(4,gint),xlbdd(4,mint),omega(3,gint)
      real*8    thts(3,gint),thtpm(3,gint),thtd(3,mint)

      save

c     Computation of rotation matrix at t=n+1 and curvatures at t=n+1
c     for Gauss points for stiffness

      ilbd = 1
      iom  = 4*gint + 4*mint + 1
      do j1 = 1,gint

c       Computation of rotation matrix and curvature at t=n+1

        do k1 = 1,3
           t1(k1) = thts(k1,j1)
           t2(k1) = thtpm(k1,j1)
        end do ! k1
        call uprm3f(xlbds(1,j1),t1,newdb(ilbd))
        call upom3f(omega(1,j1),t1,t2,newdb(iom))
        ilbd = ilbd + 4
        iom  = iom + 3
      end do ! j1

c     Computation of rotation matrix at t=n+1 for Gauss points for mass

      ilbd = 4*gint + 1
      do j1 = 1,mint
        do k1 = 1,3
          t1(k1) = thtd(k1,j1)
        end do ! k1
        call uprm3f(xlbdd(1,j1),t1,newdb(ilbd))
        ilbd = ilbd + 4
      end do ! j1

c     Transfer of new database to partition from nh2

      do j1 = 1,mct0
        hr(j1+nh2-1) = newdb(j1)
      end do ! j1

      end

      subroutine agen3f(d,ul,xl,ndf,ndm,xlbds,omega,thts,thtpm,
     &                  hn,h1,nh,isw)

c     Computes angular momentum, kinetic, potential and total energies
c     at start of new time step & prints then in corresponding output
c     files
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i     xlbds(4,gint)  : Quaternions associated with rotation matrices
c                        at Gauss points for stiffness at t=n
c i     omega(3,gint)  : Curvature at Gauss points for stiffness at t=n
c i     thts(3,gint)   : Rotation vectors at Gauss points for stiffness
c                        between t=n and t=n+alpha
c i     thtpm(3,gint)  : Derivative with respect to S of rotation
c                        vectors at Gauss points for stiffness between
c                        t=n and t=n+alpha
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: FRM3FD
c CALL     : SHP1D, UPRM3F, TD3F, PRSP3F, MSTS3F
c            PZERO, PMOVE, INT1D, INT1DL, MATQUA, MULTFD
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'
      include  'cdata.h'
      include  'ddata.h'
      include  'eldata.h'
      include  'erotas.h'
      include  'tdata.h'
      include  'ptdat6.h'

      integer   ndf,ndm, j1,k1,k2,l1,nh,isw
      real*8    d(*),ul(ndf,nen,*),xl(ndm,nel),xlbds(4,gint)
      real*8    omega(3,gint),thts(3,gint),thtpm(3,gint),hn(*),h1(*)
      real*8    shp(2,3),xjac,dx,t1(3),t2(3),fi(3),vn1(3)
      real*8    rotm1(3,3),dx1,vsnn1(3,3),wsnn1(3,3),xlbd(4)
      real*8    sg(2,5),fipr(3),omg(3),mstrn(6),mstrs(6)
      real*8    asnn1(3,3),arsnn1(3,3), con1,con2, dot
      real*8    rotnn1(3,3,3),sinn1(3,3,3),angnn1(3,3),aux1(3,3)
      real*8    rot(3),rotp(3),linm(3),angm(3),enerk,enerp,atnorm,arnorm

      save

c     Set element angular momentum and energies to zero

      do j1 = 1,3
        linm(j1) = 0.0d0
        angm(j1) = 0.0d0
      end do ! j1
      enerk  = 0.d0
      enerp  = 0.d0
      atnorm = 0.d0
      arnorm = 0.d0

c     Transfer of spatial velocities at nodes at t=n+1 and spatial
c     accelerations at nodes

      con1  = 2.d0/dt
      con2  = con1*con1
      do j1 = 1,nel
        do k1 = 1,3
          vsnn1(k1,j1)  = ul(k1,j1,4)
          asnn1(k1,j1)  = con1*vsnn1(k1,j1) - con2*ul(k1,j1,2)
          wsnn1(k1,j1)  = rvel(k1,j1,2)
          arsnn1(k1,j1) = racc(k1,j1,2)
        end do ! k1

c       Computation of rotation matrices of nodes tt t=n and t=n+1

        call quamat(xln(1,j1,3),rotnn1(1,1,j1))

c       Computation of inertia dyadic at nodes at t=n and t=n+1

        do l1 = 1,3
          do k1 = 1,3
            aux1(k1,l1) = rotnn1(k1,1,j1)*miner(1,l1)
     &                  + rotnn1(k1,2,j1)*miner(2,l1)
     &                  + rotnn1(k1,3,j1)*miner(3,l1)
          end do ! k1
        end do ! l1

        do l1 = 1,3
          do k1 = 1,3
            sinn1(k1,l1,j1) = aux1(k1,1)*rotnn1(l1,1,j1)
     &                      + aux1(k1,2)*rotnn1(l1,2,j1)
     &                      + aux1(k1,3)*rotnn1(l1,3,j1)
          end do ! k1
        end do ! l1

c       Computation of angular momentum at nodes at t=n and t=n+1

        do k1 = 1,3
          angnn1(k1,j1) = sinn1(k1,1,j1)*wsnn1(1,j1)
     &                  + sinn1(k1,2,j1)*wsnn1(2,j1)
     &                  + sinn1(k1,3,j1)*wsnn1(3,j1)
        end do ! k1

      end do ! j1

c     Computation of norm of accelerations

      do j1 = 1,nel-1
        do k1 = 1,3
          atnorm = atnorm +  asnn1(k1,j1)* asnn1(k1,j1)
          arnorm = arnorm + arsnn1(k1,j1)*arsnn1(k1,j1)
        end do ! k1
      end do ! j1

c     Gauss points for mass

      if(lobatt.eq.0) then
        call int1d (mint,sg)
      else
        call int1dl(mint,sg)
      end if

c     Loop over Gauss points for mass

      do j1 = 1,mint

c       Shape functions

        call shp1d(sg(1,j1),xl,shp,ndm,nel,xjac)
        dx  = sg(2,j1)*xjac
        dx1 = arho*dx

c       Computation of velocities at t=n+1 & current configuration &
c       its derivative with respect to S

        do k1 = 1,3
          fi(k1)  = 0.0d0
          vn1(k1) = 0.0d0
          t1(k1)  = 0.0d0
        end do ! k1
        do k1 = 1,nel
          do l1 = 1,3
            fi(l1)  = fi(l1) + shp(2,k1)*(xl(l1,k1) + ul(l1,k1,1))
            vn1(l1) = vn1(l1) + shp(2,k1)*vsnn1(l1,k1)
            t1(l1)  = t1(l1)  + shp(2,k1)*angnn1(l1,k1)
          end do ! l1
          enerk = enerk + 0.5d0*dx*shp(2,k1)
     &                  * dot(wsnn1(1,k1),angnn1(1,k1),3)
        end do ! k1

c       Computation of translational and angular momentum

        call vecp(fi,vn1,t2)
        do k1 = 1,3
          angm(k1) = angm(k1) + dx1*t2(k1) + dx*t1(k1)
          t1(k1)   = 0.0d0
        end do ! k1
        do k1 = 1,nel
          do l1 = 1,3
            t1(l1) = t1(l1) + angnn1(l1,k1)
          end do ! l1
        end do ! k1

c       Computation of linear momentum

        do k1 = 1,3
          linm(k1) = linm(k1) + dx1*vn1(k1)
        end do ! k1

c       Computation of kinetic energy

        enerk = enerk + 0.5d0*dx1*dot(vn1,vn1,3)

      end do ! j1

c     Gauss points for stiffness

      if(nint(d(182)).gt.0) then
        call int1dn(gint, sg)
      else
        call int1d(gint, sg)
      endif

c     Loop over Gauss points for stiffness

      do j1 = 1,gint

c       Shape functions

        call shp1d(sg(1,j1),xl,shp,ndm,nel,xjac)
        dx = sg(2,j1)*xjac

c       Computation of derivative with respect to S of current
c       configuration and current rotation matrix

        do k2 = 1,3
          fipr(k2) = 0.0d0
          rot(k2)  = 0.0d0
          rotp(k2) = 0.0d0
        end do ! k2

        do k1 = 1,nel
          do l1 = 1,3
            fipr(l1) = fipr(l1) + shp(1,k1)*(xl(l1,k1) + ul(l1,k1,1))
          end do ! l1
        end do ! k1

        do k1 = 1,nel
          do k2 = 1,3
            rot(k2)  = rot(k2)  + rots(k2,k1,2)*shp(2,k1)
            rotp(k2) = rotp(k2) + rots(k2,k1,2)*shp(1,k1)
          end do ! k2
        end do ! k1
        call upom3f(omega(1,j1),rot,rotp,omg)
        do k1 = 1,3
          t1(k1) =  thts(k1,j1)
          t2(k1) = thtpm(k1,j1)
        end do ! k1
        call uprm3f(xlbds(1,j1),t1,xlbd)
        call quamat(xlbd,rotm1)

c       Computation of material strains and stresses

        call msts3f(d,fipr,fipr,rotm1,rotm1,omg,omg,mstrn,mstrs,
     &              hn,h1,nh,isw)

c       Computation of potential energy

        enerp = enerp + 0.5d0*dx*dot(mstrn,mstrs,6)

      end do ! j1

c     transfer linear and angular momentum and energies

      epl(1) = epl(1) + linm(1)
      epl(2) = epl(2) + linm(2)
      epl(3) = epl(3) + linm(3)

      epl(4) = epl(4) + angm(1)
      epl(5) = epl(5) + angm(2)
      epl(6) = epl(6) + angm(3)

      epl(7) = epl(7) + enerk
      epl(8) = epl(8) + enerp

      end

      subroutine td3f(tht,ttht)

c     Computes matrix corresponding to linearized rotation vector for
c     a given rotation vector

c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i     tht(3)         : The rotation vector
c   o   ttht(3,3)      : The matrix T
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: DYNA3F
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,k
      real*8    tht(3),ttht(3,3)

c     Computes norm of increment

      ttht(1,1) =  1.d0
      ttht(1,2) = +tht(3)*0.5d0
      ttht(1,3) = -tht(2)*0.5d0
      ttht(2,1) = -tht(3)*0.5d0
      ttht(2,2) =  1.d0
      ttht(2,3) = +tht(1)*0.5d0
      ttht(3,1) = +tht(2)*0.5d0
      ttht(3,2) = -tht(1)*0.5d0
      ttht(3,3) =  1.d0

      do i = 1,3
        do k = 1,3
          ttht(i,k) = ttht(i,k) + .25d0*(tht(i)*tht(k))
        end do ! k
      end do ! i

      end

      subroutine hat3f(v,o)

c     Computes hat-matrix corresponding to a vector

c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i     v(3)         :  vector
c o   o(3,3)         : vector-hat
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: DYNA3F
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    v(3),o(3,3)

      o(1,1) =  0.0d0
      o(1,2) = -v(3)
      o(1,3) =  v(2)
      o(2,1) =  v(3)
      o(2,2) =  0.0d0
      o(2,3) = -v(1)
      o(3,1) = -v(2)
      o(3,2) =  v(1)
      o(3,3) =  0.0d0

      end

      subroutine init3f(d)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Obtains variables of common /bm3f1/,/bm3f2/ from element
c              properties
c CALLED BY : FRM3FD
c CALL :      INIR3F
c             PZERO, VECP, PMOVE
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'eqsym.h'

      integer   i,  j
      real*8    gg, d(*)

c     Some constants in commons /bm3f1/,/bm3f2/

c     if(nel.eq.2) then
        gint   = 1
        mint   = 3
c     else
c       gint   = 2
c       mint   = 4
c     endif
      if(nint(d(182)).gt.0) then      ! Nodal quadrature
        gint = nel
      endif

      if(nint(d(6)).lt.0) then
        lobatt = 1
      else
        lobatt = 0
      endif

c     Set flags

      if(neqs.eq.neq) then
        isym = 1
      else
        isym = 0
      endif

c     Set distributed loading type

      iforc  = int(d(69))

c     Computation of mct0

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

c     Material Hookean tensor

      do i = 1,6
        do j = 1,6
          mhook(j,i) = 0.0d0
        end do ! j
      end do ! i

      gg         = d(1)/(1.d0 + d(2))*0.5d0
      mhook(1,1) = d(37)*gg*d(32)
      mhook(2,2) = d(38)*gg*d(32)
      mhook(3,3) = d(1)*d(32)
      mhook(4,4) = d(1)*d(33)
      mhook(5,5) = d(1)*d(34)
      mhook(4,5) = d(1)*d(35)
      mhook(5,4) = mhook(4,5)
      mhook(6,6) = gg*d(36)
      mhook(1,6) = mhook(1,2)*d(94) - mhook(1,1)*d(95)
      mhook(6,1) = mhook(1,6)
      mhook(2,6) = mhook(2,2)*d(94) - mhook(1,2)*d(95)
      mhook(6,2) = mhook(2,6)

c     External loading

      eforce(1) = d(11)*dm
      eforce(2) = d(12)*dm
      eforce(3) = d(13)*dm

      end

      subroutine inir3f(d,xl,ndm,nel)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Initialize rotation matrices for Gauss points of
c              stiffness & mass & nodes
c CALLED BY : INIT3F
c CALL :      SHP1D
c             PZERO, VECP, PMOVE, INT1D, INT1DL
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'

      integer   ndm,nel, j1,k1,l1
      real*8    d(*),xl(ndm,nel),xlbda0(3,3),shp(2,3),xjac,xi,aa
      real*8    t1(3),t2(3),t3(3),dl,sg(2,5)

      save

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

c     Straight element. Computation of unique rotation matrix

      do j1 = 1,3
        t3(j1) = xl(j1,2) - xl(j1,1)
      end do ! j1
      dl = 1.d0/sqrt(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3))
      do j1 = 1,3
        t3(j1) = t3(j1)*dl
      end do ! j1

c     Transfer to rotqu0

      call fraftr(d,xl,ndm, t1,t2,t3)
      do k1 = 1,3
        xlbda0(k1,1) = t1(k1)
        xlbda0(k1,2) = t2(k1)
        xlbda0(k1,3) = t3(k1)
      end do ! k1

      call quatex(xlbda0,rotqu0(1,1))
      do j1 = 2,gint+mint+nel
        do k1 = 1,4
          rotqu0(k1,j1) = rotqu0(k1,1)
        end do ! k1
      end do ! j1

      return

c     Curved element

c     Gauss points for stifness

  101 if(nint(d(182)).gt.0) then
        call int1dn(gint, sg)
      else
        call int1d(gint, sg)
      endif

      do j1 = 1,gint
        call shp1d(sg(1,j1),xl,shp,ndm,nel,xjac)
        do k1 = 1,3
          t2(k1) = 0.0d0
          t3(k1) = 0.0d0
        end do ! k1

c       The third vector has direction of tangent to line of centroids

        do l1 = 1,nel
          do k1 = 1,3
            t2(k1) = t2(k1) + shp(2,l1)*xl(k1,l1)
            t3(k1) = t3(k1) + shp(1,l1)*xl(k1,l1)
          end do ! k1
        end do ! l1
        dl = 1.d0/sqrt(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3))
        do k1 = 1,3
          t3(k1) = t3(k1)*dl
        end do ! k1

c       Transfer to rotqu0

        call fraftr(d,xl,ndm,  t1,t2,t3)
        do k1 = 1,3
          xlbda0(k1,1) = t1(k1)
          xlbda0(k1,2) = t2(k1)
          xlbda0(k1,3) = t3(k1)
        end do ! k1
        call quatex(xlbda0,rotqu0(1,j1))
      end do ! j1

c     Gauss points for mass

      if(lobatt.eq.0) then
        call int1d (mint,sg)
      else
        call int1dl(mint,sg)
      end if
      do j1 = 1,mint
        call shp1d(sg(1,j1),xl,shp,ndm,nel,xjac)
        do k1 = 1,3
          t2(k1) = 0.0d0
          t3(k1) = 0.0d0
        end do ! k1

c       The third vector has direction of tangent to line of centroids

        do l1 = 1,nel
          do k1 = 1,3
            t2(k1) = t2(k1) + shp(2,l1)*xl(k1,l1)
            t3(k1) = t3(k1) + shp(1,l1)*xl(k1,l1)
          end do ! k1
        end do ! l1
        dl = 1.d0/sqrt(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3))
        do k1 = 1,3
          t3(k1) = t3(k1)*dl
        end do ! k1

c       Transfer to rotqu0

        call fraftr(d,xl,ndm,  t1,t2,t3)
        do k1 = 1,3
          xlbda0(k1,1) = t1(k1)
          xlbda0(k1,2) = t2(k1)
          xlbda0(k1,3) = t3(k1)
        end do ! k1
        call quatex(xlbda0,rotqu0(1,j1+gint))
      end do ! j1

c     Nodes

      do j1 = 1,nel
        xi = dble(j1) - 2.d0
        call shp1d(xi,xl,shp,ndm,nel,xjac)
        do k1 = 1,3
          t3(k1) = 0.0d0
        end do ! k1

c       The third vector has direction of tangent to line of centroids

        do k1 = 1,3
          t2(k1) = xl(k1,j1)
          do l1 = 1,nel
            t3(k1) = t3(k1) + shp(1,l1)*xl(k1,l1)
          end do ! l1
        end do ! k1
        dl = 1.d0/sqrt(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3))
        do k1 = 1,3
          t3(k1) = t3(k1)*dl
        end do ! k1

c       Transfer to rotqu0

        call fraftr(d,xl,ndm, t1,t2,t3)
        do k1 = 1,3
          xlbda0(k1,1) = t1(k1)
          xlbda0(k1,2) = t2(k1)
          xlbda0(k1,3) = t3(k1)
        end do ! k1
        call quatex(xlbda0,rotqu0(1,j1+gint+mint))
      end do ! j1

      end

      subroutine outp3f(d,ul,xl,s,p,ndf,ndm,nst,xlbds,omega,thts,thtpm,
     &                  hn,h1,nh,isw)

c     Output of variables
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i     xlbds(4,gint)  : Quaternions associated with rotation matrices
c                        at Gauss points for stiffness at t=n
c i     omega(3,2*int) : Curvatures at Gauss points of stiffness at t=n
c i     thts(3,gint)   : Rotation vectors at Gauss points for stiffness
c                        between t=n and t=n+alpha
c i     thtpm(3,gint)  : Derivative with respect to S of rotation
c                        vectors at Gauss points for stiffness between
c                        t=n and t=n+alpha
c-----[--.----+----.----+----.-----------------------------------------]
c CALLED BY: FRM3FD
c CALL     : SHP1D, TRCF3F, UPRM3F, UPOM3F, MSTS3F
c            PZERO, MATQUA, INT1D
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm3fd.h'
      include  'cdata.h'
      include  'ddata.h'
      include  'eldata.h'
      include  'iofile.h'

      integer   ndf,ndm,nst,nh,isw, i1,j1,k1
      real*8    d(*),ul(ndf,1),xl(ndm,1),s(nst,nst),p(nst),hn(*),h1(*)
      real*8    xlbds(4,gint),omega(3,gint),thts(3,gint),thtpm(3,gint)
      real*8    aa,sg(2,5),xnode(3,3),shp(2,3),xjac,dx,coor(3),fi(3)
      real*8    fiprim(3),rot(3),rotpm(3),xlbd(4),omg(3),rotmat(9)
      real*8    strain(6),stress(6),fp(12),btanc(6,18)

      save

c     Nodal translational displacements between n and n+1 for load incr.

      do i1 = 1,12
        fp(i1) = 0.0d0
      end do ! i1
      aa = (1.d0-theta(3))/theta(3)
      do i1 = 1,nel
        do j1 = 1,3
          xnode(j1,i1) = xl(j1,i1) + ul(j1,i1)+ul(j1,nen+i1)*aa
        end do ! j1
      end do ! i1

c     Computation & output of current configuration, strains & stresses
c     at Gauss points for stiffness at t=n+1

c     Gauss points

      if(nint(d(182)).gt.0) then
        call int1dn(gint, sg)
      else
        call int1d(gint, sg)
      endif

c     Computation and output for each Gauss point

      if(isw.eq.4) then
        write(iow,20000)
        if(ior.lt.0) then
          write(*,20000)
        endif
      endif

      do i1 = 1,gint

c       Coordinates of point

        call shp1d(sg(1,i1),xl,shp,ndm,nel,xjac)
        dx = sg(2,i1)*xjac
        do k1 = 1,3
          coor(k1) = 0.0d0
        end do ! k1
        do j1 = 1,nel
          do k1 = 1,3
            coor(k1) = coor(k1)+xl(k1,j1)*shp(2,j1)
          end do ! k1
        end do ! j1

c       Computation of displacements & derivative with respect to S

        call trcf3f(xnode,nel,shp,fi,fiprim)

c       Computation of rotation matrix and curvature at t=n+1

        do j1 = 1,3
          rot(j1)   = thts(j1,i1)
          rotpm(j1) = thtpm(j1,i1)
        end do ! j1
        call uprm3f(xlbds(1,i1),rot,xlbd)
        call upom3f(omega(1,i1),rot,rotpm,omg)

c       Computation of rotation matrix from quaternion

        call quamat(xlbd,rotmat)

c       Computation of material strains and stresses

        call msts3f(d,fiprim,fiprim,rotmat,rotmat,omg
     &             ,omg,strain,stress,hn,h1,nh,isw)

c       Output to file

        if(isw.eq.4) then

c         Output of initial coordinates, and configuration

          write(iow,20001) n,i1,coor,fi,rotmat
          if(ior.lt.0) then
            write(*,20001) n,i1,coor,fi,rotmat
          endif

c         Output of material strains and stresses and spatial curvature

          write(iow,20002) strain,stress,omg
          if(ior.lt.0) then
            write(*,20002) strain,stress,omg
          endif
        elseif(isw.eq.8) then
          call bc3f(shp,nel,fiprim,btanc)

          do j1 = 1,6
            aa = stress(j1)*dx
            do k1 = 1,12
              fp(k1) = fp(k1) + btanc(j1,k1)*aa
            end do ! k1
          end do ! kl
        endif

      end do ! i1

c     Project stress states

      if(isw.eq.8) then
        call frcn3d(fp,p,s)
      endif

20000 format(/5x,'V A R I A B L E S   A T   G A U S S   P O I N T S'/
     &        4x,50('*')/)

20001 format(5x,'Element No.',i7,': Point No.',i2/
     &       5x,'Initial Coordinates'/5x,1p,3e15.5/
     &       5x,'Current Coordinates'/5x,1p,3e15.5/
     &       5x,'Current Intrinsic Frame'/5x,3(1p,3e15.5/5x))

20002 format(5x,'Material Strains'/1p,6e13.5/
     &       5x,'Material Stresses'/1p,6e13.5/
     &       5x,'Spatial Curvature'/5x,1p,3e15.5/)

      end

      subroutine multfd(a,b,c,nra,nca,ncb)

c     Computes product of two matrices
c-----[--.----+----.----+----.-----------------------------------------]
c - ARGUMENTS
c i       a(nra,nca) - First factor of product
c i       b(nca,ncb) - Second factor of product
c   o     c(nra,ncb) - Result
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nra,nca,ncb, i,j,k
      real*8    a(nra,nca),b(nca,ncb),c(nra,ncb)

c     Computation of product

      do j = 1,nra
        do i = 1,ncb
          c(i,j) = 0.d0
          do k = 1,nca
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          end do ! k
        end do ! i
      end do ! j

      end

      subroutine fraftr(d,xl,ndm, t1,t2,t3)

      implicit   none

      include   'refnd.h'

      integer    ndm,j1
      real*8     d(*),xl(ndm,*), t1(3),t2(3),t3(3), dl,theta

c     Set transformation type

      lref    = nint(d(96))
      refx(1) = d(97)
      refx(2) = d(98)
      refx(3) = d(99)

c     Reference node

      if    (lref.eq.1) then
        do j1 = 1,3
          t2(j1) = refx(j1) - xl(j1,1)
        end do ! j1

c     Reference vector

      elseif(lref.eq.2) then
        do j1 = 1,3
          t2(j1) = refx(j1)
        end do ! j1

c     Reference polar

      elseif(lref.eq.3) then
        t2(1) = 0.5d0*(xl(1,1) + xl(1,2))
        t2(2) = 0.5d0*(xl(2,1) + xl(2,2))
        theta = atan2(t2(2),t2(1))
        dl    = sqrt(t2(1)**2 + t2(2)**2)
        t2(1) = dl*cos(theta)
        t2(2) = dl*sin(theta)
        t2(3) = 0.0d0
      endif

      call vecp(t2,t3,t1)
      dl = 1.d0/sqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
      do j1 = 1,3
         t1(j1) = t1(j1)*dl
      end do ! j1
      call vecp(t3,t1,t2)

      end
