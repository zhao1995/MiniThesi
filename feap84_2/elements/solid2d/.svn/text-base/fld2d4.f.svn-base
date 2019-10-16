c$Id:$
       subroutine fld2d4(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add radial body loading                          18/05/2007
c       2. Pass strains to stcn2z for z-zhu projections     01/01/2014
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 4-Node Energy Conserving Formulation

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ix(*)     - Global nodal connections
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         r(ndf,*)  - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'iofile.h'
      include  'prstrs.h'
      include  'strnum.h'
      include  'comblk.h'

      integer   ndf,ndm,nst,isw, i,j,i1,j1,l,lint, ix(*)
      real*8    d(*),ul(ndf,*),xl(ndm,*),s(nst,*),p(*)
      real*8    sg(3,9),shp(3,4),rho,xx1,xsj

      save

c     Data inputs

      if( isw.eq.1 ) then

        if(nint(d(16)).eq.8 .and. nint(d(19)).ne.3) ix(3) = 0

c     Tangent stiffness and residual vector

      elseif( isw.eq.3 .or. isw.eq.6 ) then

        call stif2d4(ndm,ndf,nst, d,ul,xl,ix, p,s, isw)

c     Outputs

      elseif( isw.eq.4.or.isw.eq.8.or.isw.eq.16.or.isw.eq.25) then

        call stre2d4(ndm,ndf, ix,d,ul,xl, p,s, isw)

c     Mass computation

      elseif( isw.eq.5 ) then

c       Get quadrature information

        if(nint(d(182)).gt.0) then
          call int2dn(nel,lint,sg)
        else
          l = 2
          call int2d(l,lint, sg)
        endif

        rho = d(4)

c       LOOP OVER GAUSS POINTS

        do l = 1,lint

c         Compute shape functions

          call shp2d(sg(1,l),xl,shp,xsj,ndm,4,ix,.false.)
          dm = sg(3,l)*xsj*rho

c         For each node j compute db = rho*shape*dm

          j1 = 1
          do j = 1,nel
            xx1 = shp(2,j)*dm

c           Compute a lumped mass

            p(j1)   = p(j1) + xx1
            p(j1+1) = p(j1)

c           For each node i compute mass matrix (upper triangular part)

            i1 = 1
            do i = 1,nel
              s(j1  ,i1  ) = s(j1,i1) + shp(3,i)*xx1
              s(j1+1,i1+1) = s(j1,i1)
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j

        end do ! l

c     Energy and Momenta

      elseif( isw.eq.13 ) then

        call ener2d4(ndm,ndf,nel, d,ul,xl,ix)

      end if

      end

      subroutine bmat2d4(f,shp,bb, nel)

c-----[--.----+----.----+----.-----------------------------------------]
c     Non-linear B-matrix

c     Input:
c       f(3,3)      - deformation gradient
c       shp(3,nel)  - shape functions
c       nel         - number of nodes on element
c     Output:
c       bb(4,2,nel) - strain displacement matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  a,i,j,nel
      real*8   f(3,3),shp(3,*), bb(4,2,*)

c     Loop over nodes on element

      do a = 1,nel
        do i = 1,2
          do j = 1,2

c           Normal strain terms

            bb(j,i,a) = f(i,j)*shp(j,a)
          end do ! j

c         Shear strain terms

          bb(3,i,a) = 0.0d0
          bb(4,i,a) = f(i,1)*shp(2,a) + f(i,2)*shp(1,a)
        end do ! i
      end do ! a

      end

      subroutine ener2d4(ndm,ndf,nel, d,ul,xl,ix)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Momentum/energy computation in reference coordinates

c     Input:
c       ndm           - Number of space dimensions in mesh
c       ndf           - Number of dof/node
c       nel           - Number of nodes on element
c       d(*)          - Material parameter array
c       ul(ndf,nen,*) - Nodal solution array
c       xl(ndm,*)     - Nodal coordinate array

c     Output:
c       Energy and momenta through common /ptdat6/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'cdata.h'
      include 'eltran.h'
      include 'ptdat6.h'
      include 'tdata.h'

      integer  a, i,j,l, lint, ndm,ndf,nel, ix(4)
      real*8   detf, xsj, dvol, con1,con2, kineng, poteng
      real*8   fa(3,3),f0(3,3), shp(3,4,4), ss(6), dm(6,6)
      real*8   ee(6),  sg(3,4), vel(2), mom(2,4)
      real*8   d(*), ul(ndf,nen,*), xl(ndm,*)

      save

c     Initialize

      do a = 1,nel
        do j = 1,2
          mom(j,a)  = 0.0d0
        end do ! j
      end do ! a

      kineng = 0.0d0
      poteng = 0.0d0

c     Get quadrature information

      if(nint(d(182)).gt.0) then
        call int2dn(nel,lint,sg)
      else
        l = 2
        call int2d(l,lint, sg)
      endif

c     Loop over quadrature points

      do l = 1,lint

        call shp2d(sg(1,l),xl,shp(1,1,l),xsj,ndm,4,ix,.false.)

        dvol = xsj*sg(3,l)

c       Compute deformation gradient and increment

        call fmat2d4(shp(1,1,l),ul,ndf,nel,nen, fa,f0,detf)

c       Compute Green-Lagrange strain

        call emat3d4(1.d0, fa,f0, ee)

        call modl3d4(d,detf,fa, dm,ss,ee,.false.)

c       Accumlate Potential energy for element

        poteng = poteng + (ss(1)*ee(1) + ss(2)*ee(2) + ss(3)*ee(3)
     &                  +  ss(4)*ee(4))*dvol

c       Form velocity at gauss point

        con2 = d(4) * dvol
        con1 = con1 * d(7)  ! Constant factor

c       Inertial contribution from consistent mass approximation

        do i = 1,2
          vel(i) = 0.0d0
          do a = 1,nel
            vel(i) = vel(i) + shp(3,a,l)*ul(i,a,4)
          end do ! a

          do a = 1,nel
            mom(i,a) = mom(i,a) + shp(3,a,l)*vel(i)*con1
          end do ! a

c         Accumulate kinetic energy for element

          kineng = kineng + vel(i)*vel(i)*con1

        end do ! i

c       Inertial contribution from lumped mass approximation

        con2 = con2 - con1       ! Lumped     factor
        do a = 1,nel
          con1 = con2*shp(3,a,l) ! Lumped times shape function
          do i = 1,2
            mom(i,a) = mom(i,a) + ul(i,a,4)*con1
            kineng   = kineng   + ul(i,a,4)*ul(i,a,4)*con1
          end do ! i
        end do ! a

      end do ! l

c     Accumulate total momenta (xsm = deformed x)

      do a = 1,nel

c       Linear momentum

        epl(1) = epl(1) + mom(1,a)
        epl(2) = epl(2) + mom(2,a)

c       Angular momentum

        epl(6) = epl(6) + (xl(1,a) + ul(1,a,1))*mom(2,a)
     &                  - (xl(2,a) + ul(2,a,1))*mom(1,a)
      end do ! a

c     Accumulate total kinetic and stored energy

      epl(7) = epl(7) + 0.5d0*kineng
      epl(8) = epl(8) + 0.5d0*poteng

      end

      subroutine fmat2d4(shp,ul,ndf,nel,nen, f,df,detf)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Form Deformation Gradient

c     Inputs:

c       shp(3,4)      - shape function derivatives: N_I,i = shp(i,I)
c       ul(ndf,nen,1) - nodal displacements
c       ul(ndf,nen,2) - increment to nodal displacements
c       ndf           - number of degree of freedoms at node
c       nel           - number of nodes on element
c       nen           - max number of nodes on element

c     Outputs:

c       f(3,3)      - deformation gradient
c       df(3,3)     - incremental deformation gradient
c       detf        - determinant of deformation gradient
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  ndf,nel,nen, i,j,k
      real*8   shp(3,4),ul(ndf,nen,*), f(3,3),df(3,3), detf

      do i = 1,3
        do j = 1,3
          f(i,j)  = 0.0d0
          df(i,j) = 0.0d0
        end do ! j
        f(i,i) = 1.0d0
      end do ! i
      do i = 1,2
        do j = 1,2
          do k = 1,nel
            f(i,j)  = f(i,j)  + ul(i,k,1)*shp(j,k)
            df(i,j) = df(i,j) + ul(i,k,2)*shp(j,k)
          end do ! k
        end do ! j
      end do ! i

      detf = f(1,1)*f(2,2) - f(1,2)*f(2,1)

      end

      subroutine gstf2d4(ss,shp,rhocdv,rholdv,ctan,ndf,nel,nst, s)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Geometric stiffness and consistent mass contribution

c     Input:
c       ss(6)      - Stress times d_volume
c       shp(3,4)   - Shape function derivatives
c       ctan(3)    - Algorithm stiffness parameters
c       rhocdv     - Mass density times volume - consistent
c       rholdv     - Mass density times volume - lumped
c       ndf        - Number of degree-of-freedoms per node
c       nel        - Number of nodes on element
c       nst        - Size of element stiffness

c     Output:
c       s(nst,nst) - geometric stiffness
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  ndf,nel,nst, a,b, ia, ib, i
      real*8   ss(6),shp(3,4),s(nst,nst),t(2),ctan(3),kgeom
      real*8   rholdv,rhocdv,dmas

      ib = 0
      do b = 1,nel

c       Lumped mass contribution

        dmas = rholdv*shp(3,b)*ctan(3)
        do i = 1,2
          s(ib+i,ib+i) = s(ib+i,ib+i) + dmas
        end do ! i
        t(1) = (ss(1)*shp(1,b)+ss(4)*shp(2,b))*ctan(1)
        t(2) = (ss(4)*shp(1,b)+ss(2)*shp(2,b))*ctan(1)

c       Consistent mass contribution

        dmas = rhocdv*shp(3,b)*ctan(3)

        ia = 0
        do a = 1,nel
          kgeom = shp(1,a)*t(1) + shp(2,a)*t(2) + shp(3,a)*dmas
          do i = 1,2
            s(ia+i,ib+i) = s(ia+i,ib+i) + kgeom
          end do ! i
          ia = ia + ndf
        end do ! a
        ib = ib + ndf
      end do ! b

      end

      subroutine stif2d4(ndm,ndf,nst, d,ul,xl,ix, p,s, isw)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Stiffness/residual formulation in reference coordinates

c     Input:
c       ndm           - Number of coordinates in mesh per node
c       ndf           - Number of degree-of-freedoms per node
c       nst           - Size of element stiffness
c       d(*)          - Material parameter array
c       ul(ndf,nen,*) - Nodal solution array
c       xl(ndm,*)     - Nodal coordinate array
c       ix(*)         - Global node numbers on element
c       isw           - Switch to control computations

c     Output:
c       p(ndf,*)      - Residual
c       s(nst,nst)    - Stiffness
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'cdata.h'
      include 'eldata.h'
      include 'elplot.h'
      include 'eltran.h'
      include 'part0.h'
      include 'pmod2d.h'
      include 'ptdat6.h'
      include 'rdata.h'
      include 'tdata.h'

      integer  a,b, i,j,k,l, ia,jb, isw, lint, ndm,ndf,nst, ix(4)
      real*8   fa(3,3),shp(3,4,4), bb(4,2,4), bd(2,4), ss(6), ds(6,6)
      real*8   dvol(4),ee(6),sg(3,4),ace(2),body(3),bdy(3), detf, xsj
      real*8   rhol,rhoc, con1,con2
      real*8   df(3,3),f0(3,3),f1(3,3),bba(4,2,4),xr(2)
      real*8   d(*), ul(ndf,nen,*), xl(ndm,*), p(ndf,*), s(nst,nst)

      save

c     Compute body forces

      call sbodyf(d, body)
      bdy(1) = body(1)
      bdy(2) = body(2)

      if((d(7).ge.0.0 .or. d(183).ne.0.0d0) .and.
     &           (ndfo(1).gt.0 .or. shflg)) then
        rhol = d(4)
        rhoc = rhol * d(7)  ! Consistent
        rhol = rhol - rhoc  ! Lumped
      else
        rhol = 0.0d0
        rhoc = 0.0d0
      endif

c     Get quadrature information

      if(nint(d(182)).gt.0) then
        call int2dn(nel,lint,sg)
      else
        l = 2
        call int2d(l,lint, sg)
      endif

c     Loop over quadrature points

      do l = 1,lint

        call shp2d(sg(1,l),xl,shp(1,1,l),xsj,ndm,4,ix,.false.)
        dvol(l) = xsj*sg(3,l)

c       Compute deformation gradient and increment

        call fmat2d4(shp(1,1,l),ul,ndf,nel,nen, fa,df,detf)

        con1 = (1.d0 - ctan(1))/ctan(1)
        do j = 1,3
          do i = 1,3
            f1(i,j) = fa(i,j) + con1*df(i,j)
            f0(i,j) = fa(i,j) -      df(i,j)
          end do ! i
        end do ! j

c       Compute strain-displacement matrices

        call bmat2d4(fa,shp(1,1,l),bba,nel)
        call bmat2d4(f1,shp(1,1,l),bb ,nel)

c       Compute Green-Lagrange strain

        call emat3d4(ctan(1), f1,f0, ee)

c       Compute stresses and moduli

        call modl3d4(d,detf,f1, ds,ss,ee,.true.)

c       Save stress/strain for tplot use

        i = 12*(l-1)
        do j = 1,6
          tt(j+i  ) = ss(j)
          tt(j+i+6) = ee(j)
        end do ! j

c       Multiply stress by volume element

        con1 = dvol(l)*ctan(1)
        do i = 1,6
          ss(i) = ss(i)*dvol(l)
          do j = 1,6
            ds(i,j) = ds(i,j)*con1
          end do ! j
        end do ! i

c       Form acceleration and body loading at gauss point

        if(nint(d(69)).eq.5) then
          do i = 1,2
            xr(i) = 0.0d0
            do a = 1,nel
              xr(i) = xr(i) + shp(3,a,l)*xl(i,a)
            end do ! a
          end do ! i
          con1 = sqrt(xr(1)**2 + xr(2)**2)
          if(con1.gt.0.0d0) then
            bdy(1) = (xr(1)*body(1) - xr(2)*body(2))/con1
            bdy(2) = (xr(2)*body(1) + xr(1)*body(2))/con1
          endif
        endif

        con1 = rhoc*dvol(l)
        do i = 1,2
          ace(i) = 0.0d0
          do a = 1,nel
            ace(i) = ace(i) + shp(3,a,l)*ul(i,a,5)
          end do ! a
          ace(i) = ace(i)*con1 - bdy(i)*dvol(l)
        end do ! i

c       Compute residual

        con2 = rhol*dvol(l)
        do a = 1,nel
          do i = 1,2
            do j = 1,4
              p(i,a) = p(i,a) - bba(j,i,a)*ss(j)
            end do ! j
            p(i,a) = p(i,a) - shp(3,a,l)*(ace(i) - con2*ul(i,a,5))
          end do ! i
        end do ! a

c       Compute tangent

        if(isw.eq.3) then

c         Form geometric and inertial parts of tangent stiffness

          if(gflag) then
            call gstf2d4(ss,shp(1,1,l),con1,con2,ctan,ndf,nel,nst, s)
          endif

c         Form material part of tangent

          ia = 0
          do a = 1,nel

            do i = 1,2
              do j = 1,4
                bd(i,j) = 0.0d0
                do k = 1,4
                  bd(i,j) = bd(i,j) + bba(k,i,a)*ds(k,j)
                end do ! k
              end do ! j
            end do ! i

            jb = 0
            do b = 1,nel

              do i = 1,2
                do j = 1,2
                  do k = 1,4
                    s(ia+i,jb+j) = s(ia+i,jb+j) + bd(i,k)*bb(k,j,b)
                  end do ! k
                end do ! j
              end do ! i

              jb = jb + ndf
            end do ! b

            ia = ia + ndf
          end do ! a
        end if

      end do ! l

c     Multiply by thickness if not unity

      if((isw.eq.3 .or.isw.eq.6) .and. d(14).ne.1.d0) then

        do j = 1,nst
          do i = 1,nst
            s(i,j) = s(i,j)*d(14)
          end do ! i
        end do ! j
        do j = 1,nel
          do i = 1,ndf
            p(i,j) = p(i,j)*d(14)
          end do ! i
        end do ! j

      endif

      end

      subroutine stre2d4(ndm,ndf, ix, d,ul,xl, dt,st, isw)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Stress output in reference coordinates

c     Input:
c       ndm           - Number of coordinates in mesh per node
c       ndf           - Number of degree-of-freedoms per node
c       ix(*)         - Nodal connection list for element
c       d(*)          - Material parameter array
c       ul(ndf,nen,*) - Nodal solution array
c       xl(ndm,*)     - Nodal coordinate array
c       isw           - Switch to control computations

c     Output:
c       dt(*)         - Projection weight vector
c       st(nen,*)     - Projection array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'bdata.h'
      include 'cdata.h'
      include 'eldata.h'
      include 'iofile.h'
      include 'strnum.h'

      integer  i,j,k,l, isw, lint, ndm,ndf, ix(*)

      real*8   dt(nen), st(nen,*), d(*), ul(ndf,*), xl(ndm,*)
      real*8   ee0(6), ee(6,4), ss0(6), ss(9,4), dd(6,6), weng(4)
      real*8   f(3,3,2,4), finv(3,3,4), detf(2,4), df(3,3), spr(3)
      real*8   shp(3,9,4), dvol(4), sg(3,4),xsj, xx1,xx2

      save

      data     finv / 36*0.0d0 /

c     Get quadrature information and initialize

      if(nint(d(182)).gt.0) then
        call int2dn(nel,lint,sg)
      else
        l = 2
        call int2d(l,lint, sg)
      endif

      do i = 1,6
        ss0(i) = 0.0d0
        ee0(i) = 0.0d0
      end do ! i

      xx1 = 0.0d0
      xx2 = 0.0d0

c     Loop over quadrature points

      do l = 1,lint

        call shp2d(sg(1,l),xl,shp(1,1,l),xsj,ndm,4,ix,.false.)
        dvol(l) = xsj*sg(3,l)

c       Compute deformation gradient and increment

        call fmat2d4(shp(1,1,l),ul,ndf,nel,nen,
     &               f(1,1,1,l),df,detf(1,l))

c       Compute Green-Lagrange strain

        call emat3d4(1.d0, f,df, ee)

        call modl3d4(d,detf(1,l),f(1,1,1,l), dd,ss(1,l),ee(1,l),.false.)

        if(isw.eq.4) then

          do i = 1,4
            ss0(i) = ss0(i) + 0.25d0*ss(i,l)
            ee0(i) = ee0(i) + 0.25d0*ee(i,l)
          end do ! i

c         Compute coordinates for outputs

          do i=1,nel
            xsj = 0.25d0*shp(3,i,l)
            xx1 = xx1 + xsj*xl(1,i)
            xx2 = xx2 + xsj*xl(2,i)
          end do ! i

c       Store energy for material force computation

        elseif(isw.eq.16) then

          weng(l) = ss(1,l)*ee(1,l) + ss(2,l)*ee(2,l) + ss(3,l)*ee(3,l)
     &            + ss(4,l)*ee(4,l)

c         Compute inverse deformation gradient

          finv(1,1,l) =  f(2,2,1,l)/detf(1,l)
          finv(1,2,l) = -f(1,2,1,l)/detf(1,l)
          finv(2,1,l) = -f(2,1,1,l)/detf(1,l)
          finv(2,2,l) =  f(1,1,1,l)/detf(1,l)
          finv(3,3,l) =  f(3,3,1,l)/detf(1,l)

        end if

      end do ! l

c     Output stresses

      if(isw.eq.4) then

        mct = mct - 3
        if(mct.le.0) then
          write(iow,2001) o,head
          if(ior.lt.0) write(*,2001) o,head
          mct = 51
        endif

c       Output quadrature values

        if(qoutfl) then

          do l = 1,lint

            xx1 = 0.0d0
            xx2 = 0.0d0
            do i=1,nel
              xx1 = xx1 + shp(3,i,l)*xl(1,i)
              xx2 = xx2 + shp(3,i,l)*xl(2,i)
            end do ! i

            call pstr2d(ss(1,l),spr)

            write(iow,2002) n,ma,(ss(j,l),j=1,4),(ee(j,l),j=1,4),
     &                      xx1,xx2,spr
            if(ior.lt.0) then
              write(*,2002) n,ma,(ss(j,l),j=1,4),(ee(j,l),j=1,4),
     &                      xx1,xx2,spr
            end if
          end do ! l

c       Output averaged values

        else

          call pstr2d(ss0,spr)

          write(iow,2002) n,ma,(ss0(j),j=1,4),(ee0(j),j=1,4),xx1,xx2,spr
          if(ior.lt.0) then
            write(*,2002) n,ma,(ss0(j),j=1,4),(ee0(j),j=1,4),xx1,xx2,spr
          end if

        endif

c     Project stresses

      elseif(isw.eq.8) then

c       Perform integration

        do l = 1,lint

          do i = 1,nel
            xsj   = shp(3,i,l)*dvol(l)
            dt(i) = dt(i) + xsj
            do k = 1,4
              st(i,k  ) = st(i,k  ) + ss(k,l)*xsj
              st(i,k+6) = st(i,k+6) + ee(k,l)*xsj
            end do ! k
          end do ! i

        end do ! l

        iste = 12

c     Compute fracture indices

      elseif(isw.eq.16) then

        call pfrac2f(f,detf,ss,weng, shp,dvol, dt, lint,ndf,ndm,2)

c     Compute Z-Z projections

      elseif(isw.eq.25) then

        call stcn2z(xl,ss,ee,shp,dvol,lint,ndm,nel,9)

      end if

2001  format(a1,20a4//5x,'Element Stresses'//
     & '  Elmt  Matl   11-stress   22-stress   33-stress   12-stress'/
     & '               11-strain   22-strain   33-strain   12-strain'/
     & '     1-coord     2-coord    1-stress    2-stress       angle')

2002  format(2i6,1p,4e12.3/12x,1p,4e12.3/1p,5e12.3)

      end
