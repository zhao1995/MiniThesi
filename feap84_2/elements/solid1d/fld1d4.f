c$Id:$
      subroutine fld1d4(d,ul,xl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set quadrature for d(5).eq.0                     26/03/2009
c       2. Use 'quard1d' and 'interp1d' for solution        01/04/2009
c       3. Add flag to call list                            03/03/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 2-Node Energy Conserving Formulation

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
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
      include  'qudshp.h'
      include  'strnum.h'
      include  'comblk.h'

      integer   ndf,ndm,nst,isw, i,j,i1,j1,l
      real*8    rho,xx1
      real*8    d(*),ul(ndf,*),xl(ndm,*),s(nst,*),p(*)

      save

c     Tangent stiffness and residual vector

      if( isw.eq.3 .or. isw.eq.6 ) then

        call stif1d4(ndm,ndf,nst, d,ul,xl, p,s, isw)

c     Outputs

      elseif( isw.eq.4.or.isw.eq.8.or.isw.eq.16.or.isw.eq.25) then

        call stre1d4(ndm,ndf,nst, d,ul,xl, p,s, isw)

c     Mass computation

      elseif( isw.eq.5 ) then

c       Get quadrature information

        call quadr1d(d)

        rho = d(4)

c       LOOP OVER GAUSS POINTS

        do l = 1,lint

c         Compute shape functions

          call interp1d(l, xl, ndm,nel,.false.)
          dm = jac(l)*rho

c         For each node j compute db = rho*shape*dm

          j1 = 1
          do j = 1,nel
            xx1 = shp1(2,j,l)*dm

c           Compute a lumped mass

            p(j1)   = p(j1) + xx1

c           For each node i compute mass matrix (upper triangular part)

            i1 = 1
            do i = 1,nel
              s(j1  ,i1  ) = s(j1,i1) + shp1(2,i,l)*xx1
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j

        end do ! l

c     Energy and Momenta

      elseif( isw.eq.13 ) then

        call ener1d4(ndm,ndf,nel, d,ul,xl)

      end if

      end

      subroutine bmat1d4(f,shp1,bb, nel)

c-----[--.----+----.----+----.-----------------------------------------]
c     Non-linear B-matrix

c     Input:
c       f(3,3)      - deformation gradient
c       shp1(2,nel) - shape functions
c       nel         - number of nodes on element
c     Output:
c       bb(4,nel) - strain displacement matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  a,nel
      real*8   f(3,3),shp1(2,*), bb(*)

c     Loop over nodes on element

      do a = 1,nel
        bb(a) = f(1,1)*shp1(1,a)
      end do ! a

      end

      subroutine ener1d4(ndm,ndf,nel, d,ul,xl)

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
      include 'qudshp.h'
      include 'tdata.h'

      integer  a, l, ndm,ndf,nel
      real*8   detf, con1,con2, kineng, poteng
      real*8   fa(3,3),f0(3,3), ss(6), dm(6,6)
      real*8   ee(6),  vel, mom(4)
      real*8   d(*), ul(ndf,nen,*), xl(ndm,*)

      save

c     Initialize

      do a = 1,nel
        mom(a)  = 0.0d0
      end do ! a

      kineng = 0.0d0
      poteng = 0.0d0

c     Get quadrature information

      call quadr1d(d)

c     Loop over quadrature points

      do l = 1,lint

        call interp1d(l, xl, ndm,nel,.false.)

c       Compute deformation gradient and increment

        call fmat1d4(shp1(1,1,l),ul,ndf,nel,nen, fa,f0,detf)

c       Compute Green-Lagrange strain

        call emat3d4(1.d0, fa,f0, ee)

        call modl3d4(d,detf,fa, dm,ss,ee,.false.)

c       Accumlate Potential energy for element

        poteng = poteng + ss(1)*ee(1)*jac(l)

c       Form velocity at gauss point

        con2 = d(4) * jac(l)
        con1 = con1 * d(7)  ! Constant factor

c       Inertial contribution from consistent mass approximation

        vel = 0.0d0
        do a = 1,nel
          vel = vel + shp1(2,a,l)*ul(1,a,4)
        end do ! a

        do a = 1,nel
          mom(a) = mom(a) + shp1(2,a,l)*vel*con1
        end do ! a

c       Accumulate kinetic energy for element

        kineng = kineng + vel*vel*con1

c       Inertial contribution from lumped mass approximation

        con2 = con2 - con1       ! Lumped     factor
        do a = 1,nel
          con1 = con2*shp1(2,a,l) ! Lumped times shape function
          mom(a) = mom(a) + ul(1,a,4)*con1
          kineng = kineng + ul(1,a,4)*ul(1,a,4)*con1
        end do ! a

      end do ! l

c     Accumulate total momenta (xsm = deformed x)

      do a = 1,nel
        epl(1) = epl(1) + mom(a)
      end do ! a

c     Accumulate total kinetic and stored energy

      epl(7) = epl(7) + 0.5d0*kineng
      epl(8) = epl(8) + 0.5d0*poteng

      end

      subroutine fmat1d4(shp1,ul,ndf,nel,nen, f,df,detf)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Form Deformation Gradient

c     Inputs:

c       shp1(2,*)      - shape function derivatives: N_I,i = shp1(i,I)
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

      integer  ndf,nel,nen, i,j
      real*8   shp1(2,*),ul(ndf,nen,*), f(3,3),df(3,3), detf

      do i = 1,3
        do j = 1,3
          f(i,j)  = 0.0d0
          df(i,j) = 0.0d0
        end do ! j
        f(i,i) = 1.0d0
      end do ! i
      do j = 1,nel
        f(1,1)  = f(1,1)  + ul(1,j,1)*shp1(1,j)
        df(1,1) = df(1,1) + ul(1,j,2)*shp1(1,j)
      end do ! j

      detf = f(1,1)

      end

      subroutine gstf1d4(ss,shp1,rhocdv,rholdv,ctan,ndf,nel,nst, s)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Geometric stiffness and consistent mass contribution

c     Input:
c       ss(6)      - Stress times d_volume
c       shp1(2,4)   - Shape function derivatives
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

      integer  ndf,nel,nst, a,b, ia, ib
      real*8   ss(6),shp1(2,*),s(nst,nst),t(2),ctan(3)
      real*8   rholdv,rhocdv,dmas

      ib = 1
      do b = 1,nel

c       Lumped mass contribution

        dmas = rholdv*shp1(2,b)*ctan(3)
        s(ib,ib) = s(ib,ib) + dmas
        t(1) = ss(1)*shp1(1,b)*ctan(1)

c       Consistent mass contribution

        dmas = rhocdv*shp1(2,b)*ctan(3)

        ia = 1
        do a = 1,nel
          s(ia,ib) = s(ia,ib) + shp1(1,a)*t(1) + shp1(2,a)*dmas
          ia = ia + ndf
        end do ! a
        ib = ib + ndf
      end do ! b

      end

      subroutine stif1d4(ndm,ndf,nst, d,ul,xl, p,s, isw)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Stiffness/residual formulation in reference coordinates

c     Input:
c       ndm           - Number of coordinates in mesh per node
c       ndf           - Number of degree-of-freedoms per node
c       nst           - Size of element stiffness
c       d(*)          - Material parameter array
c       ul(ndf,nen,*) - Nodal solution array
c       xl(ndm,*)     - Nodal coordinate array
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
      include 'qudshp.h'
      include 'rdata.h'
      include 'tdata.h'

      integer  a,b, i,j,k,l, ia,jb, isw, ndm,ndf,nst
      real*8   fa(3,3), bb(4,4), bd(4), ss(6), ds(6,6)
      real*8   ee(6),ace, detf, rhol,rhoc, con1,con2
      real*8   df(3,3),f0(3,3),f1(3,3),bba(4,4), body(3)
      real*8   d(*), ul(ndf,nen,*), xl(ndm,*), p(ndf,*), s(nst,nst)

      save

      call sbodyf(d, body)

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

      call quadr1d(d)

c     Loop over quadrature points

      do l = 1,lint

        call interp1d(l, xl, ndm,nel,.false.)

c       Compute deformation gradient and increment

        call fmat1d4(shp1(1,1,l),ul,ndf,nel,nen, fa,df,detf)

        con1 = (1.d0 - ctan(1))/ctan(1)
        do j = 1,3
          do i = 1,3
            f1(i,j) = fa(i,j) + con1*df(i,j)
            f0(i,j) = fa(i,j) -      df(i,j)
          end do ! i
        end do ! j

c       Compute strain-displacement matrices

        call bmat1d4(fa,shp1(1,1,l),bba,nel)
        call bmat1d4(f1,shp1(1,1,l),bb ,nel)

c       Compute Green-Lagrange strain

        call emat3d4(ctan(1), f1,f0, ee)

c       Compute stresses and moduli

        call modl3d4(d,detf,f1, ds,ss,ee,.true.)

c       Save stress/strain for tplot use

        i = 12*(l-1)
        do j = 1,3
          tt(j+i  ) = ss(j)
          tt(j+i+6) = ee(j)
        end do ! j

c       Multiply stress by volume element

        con1 = jac(l)*ctan(1)
        do i = 1,3
          ss(i) = ss(i)*jac(l)
          do j = 1,3
            ds(i,j) = ds(i,j)*con1
          end do ! j
        end do ! i

c       Form acceleration at gauss point

        con1 = rhoc*jac(l)
        ace = 0.0d0
        do a = 1,nel
          ace = ace + shp1(2,a,l)*ul(1,a,5)
        end do ! a
        ace = ace*con1 - body(1)*jac(l)

c       Compute residual

        con2 = rhol*jac(l)
        do a = 1,nel
          do j = 1,4
            p(1,a) = p(1,a) - bba(j,a)*ss(j)
          end do ! j
          p(1,a) = p(1,a) - shp1(2,a,l)*(ace - con2*ul(1,a,5))
        end do ! a

c       Compute tangent

        if(isw.eq.3) then

c         Form geometric and inertial parts of tangent stiffness

          if(gflag) then
            call gstf1d4(ss,shp1(1,1,l),con1,con2,ctan,ndf,nel,nst, s)
          endif

c         Form material part of tangent

          ia = 1
          do a = 1,nel

            do j = 1,4
              bd(j) = 0.0d0
              do k = 1,4
                bd(j) = bd(j) + bba(k,a)*ds(k,j)
              end do ! k
            end do ! j

            jb = 1
            do b = 1,nel
              do k = 1,4
                s(ia,jb) = s(ia,jb) + bd(k)*bb(k,b)
              end do ! k
              jb = jb + ndf
            end do ! b
            ia = ia + ndf
          end do ! a
        end if

      end do ! l

      end

      subroutine stre1d4(ndm,ndf,nnp, d,ul,xl, dt,st, isw)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Stress output in reference coordinates

c     Input:
c       ndm           - Number of coordinates in mesh per node
c       ndf           - Number of degree-of-freedoms per node
c       nnp           - Dimension of stress array
c       d(*)          - Material parameter array
c       ul(ndf,nen,*) - Nodal solution array
c       xl(ndm,*)     - Nodal coordinate array
c       isw           - Switch to control computations

c     Output:
c       dt(*)         - Projection weight vector
c       st(nnp,*)     - Projection array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'bdata.h'
      include 'cdata.h'
      include 'eldata.h'
      include 'iofile.h'
      include 'qudshp.h'

      integer  i,j,k,l, isw, ndm,ndf,nnp

      real*8   dt(nnp), st(nnp,*), d(*), ul(ndf,*), xl(ndm,*)
      real*8   ee0(6), ee(6), ss0(6), ss(9,4), dd(6,6), weng(4)
      real*8   f(3,3,2,4), finv(3,3,4), detf(2,4), df(3,3)
      real*8   xsj,xx1

      save

      data     finv / 36*0.0d0 /

c     Get quadrature information and initialize

      call quadr1d(d)

      do i = 1,6
        ss0(i) = 0.0d0
        ee0(i) = 0.0d0
      end do ! i

      xx1 = 0.0d0

c     Loop over quadrature points

      do l = 1,lint

        call interp1d(l, xl, ndm,nel,.false.)

c       Compute deformation gradient and increment

        call fmat1d4(shp1(1,1,l),ul,ndf,nel,nen,
     &               f(1,1,1,l),df,detf(1,l))

c       Compute Green-Lagrange strain

        call emat3d4(1.d0, f,df, ee)

        call modl3d4(d,detf(1,l),f(1,1,1,l), dd,ss(1,l),ee,.false.)

        if(isw.eq.4) then

          do i = 1,3
            ss0(i) = ss0(i) + 0.5d0*ss(i,l)
            ee0(i) = ee0(i) + 0.5d0*ee(i)
          end do ! i

c         Compute coordinates for outputs

          do i=1,nel
            xx1 = xx1 + 0.5d0*shp1(2,i,l)*xl(1,i)
          end do ! i

c       Store energy for material force computation

        elseif(isw.eq.16) then

          weng(l) = ss(1,l)*ee(1)

c         Compute inverse deformation gradient

          finv(1,1,l) =  1.d0/f(1,1,1,l)
          finv(2,2,l) =  1.d0
          finv(3,3,l) =  1.d0/f(3,3,1,l)

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

c       Compute principal stresses

        write(iow,2002) n,ma,xx1,(ss0(j),j=1,3),(ee0(j),j=1,3)
        if(ior.lt.0) then
          write(*,2002) n,ma,xx1,(ss0(j),j=1,3),(ee0(j),j=1,3)
        end if

c     Project stresses

      elseif(isw.eq.8) then

c       Perform integration

        do l = 1,lint

          do i = 1,nel
            xsj   = shp1(2,i,l)*jac(l)
            dt(i) = dt(i) + xsj
            do k = 1,3
              st(i,k) = st(i,k) + ss(k,l)*xsj
            end do ! k
          end do ! i

        end do ! l

c     Compute fracture indices

      elseif(isw.eq.16) then

        call pfrac1f(f,detf,ss,weng, shp1,jac, dt, lint,ndf,ndm,2)

c     Compute Z-Z projections

      elseif(isw.eq.25) then

        call stcn1z(xl,ss,shp1,jac,lint,ndm,nel,9)

      end if

2001  format(a1,20a4//5x,'Element Stresses'//
     & '  Elmt  Matl   11-stress   22-stress   33-stress   12-stress'/
     & '               11-strain   22-strain   33-strain   12-strain'/
     & '     1-coord     2-coord    1-stress    2-stress       angle')

2002  format(2i6,1p,4e12.3/12x,1p,4e12.3/1p,5e12.3)

      end
