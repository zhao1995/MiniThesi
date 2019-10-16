c$Id:$
      subroutine fld3d4(d,ul,xl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: 8-Node Energy Conserving Formulation

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
c         p(ndf,*)  - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'iofile.h'
      include  'prstrs.h'
      include  'strnum.h'
      include  'comblk.h'

      integer   ndf,ndm,nst,isw
      real*8    d(*),ul(ndf,*),xl(ndm,*),s(nst,*),p(*)

      save

c     Tangent stiffness and residual vector

      if( isw.eq.3 .or. isw.eq.6 ) then

        call stif3d4(ndm,ndf,nst, d,ul,xl, p,s, isw)

c     Outputs

      elseif( isw.eq.4 .or. isw.eq.8 ) then

        call stre3d4(ndm,ndf, d,ul,xl, p,s, isw)

c     Energy and Momenta

      elseif( isw.eq.13 ) then

        call ener3d4(ndm,ndf,nel, d,ul,xl)

      end if

      end

      subroutine bmat3d4(f,shp,bb, nel)

c-----[--.----+----.----+----.-----------------------------------------]
c     Non-linear B-matrix

c     Input:
c       f(3,3) - deformation gradient
c       shp(4,nel) -shape functions
c       nel        - number of nodes on element
c     Output:
c       bb(6,3,nel) - strain displacement matrix
c-----[--.----+----.----+----.-----------------------------------------]

      implicit none

      integer  a,i,j,k,nel
      real*8   f(3,3),shp(4,*), bb(6,3,*)

c     Loop over nodes on element

      do a = 1,nel

        do i = 1,3
          k = mod(i,3) + 1
          do j = 1,3

c           Normal strain terms

            bb(i  ,j,a) = f(j,i)*shp(i,a)

c           Shear strain terms

            bb(i+3,j,a) = f(j,i)*shp(k,a) + f(j,k)*shp(i,a)

          end do ! j
        end do ! i
      end do ! a

      end

      subroutine emat3d4(al1, f,f0, ee)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Green-Lagrange strain computation

c     Input:
c       al1     - alpha for interpolation
c       f(3,3)  - deformation gradient at t+n+1
c       f0(3,3) - deformation gradient at t+n
c     Output:
c       ee(6)   - Green-Lagrange strain in matrix form
c-----[--.----+----.----+----.-----------------------------------------]

      implicit none

      integer  i
      real*8   f(*),f0(*), ee(*), ee1(6), ee0(6), al1,al0

c     Compute Green-Lagrange strains at t_n+1

      ee1(1) = 0.5d0*(f(1)*f(1) + f(2)*f(2) + f(3)*f(3)) - 0.5d0
      ee1(2) = 0.5d0*(f(4)*f(4) + f(5)*f(5) + f(6)*f(6)) - 0.5d0
      ee1(3) = 0.5d0*(f(7)*f(7) + f(8)*f(8) + f(9)*f(9)) - 0.5d0
      ee1(4) =        f(1)*f(4) + f(2)*f(5) + f(3)*f(6)
      ee1(5) =        f(4)*f(7) + f(5)*f(8) + f(6)*f(9)
      ee1(6) =        f(7)*f(1) + f(8)*f(2) + f(9)*f(3)

c     Compute Green-Lagrange strains at t_n

      ee0(1) = 0.5d0*(f0(1)*f0(1) + f0(2)*f0(2) + f0(3)*f0(3)) - 0.5d0
      ee0(2) = 0.5d0*(f0(4)*f0(4) + f0(5)*f0(5) + f0(6)*f0(6)) - 0.5d0
      ee0(3) = 0.5d0*(f0(7)*f0(7) + f0(8)*f0(8) + f0(9)*f0(9)) - 0.5d0
      ee0(4) =        f0(1)*f0(4) + f0(2)*f0(5) + f0(3)*f0(6)
      ee0(5) =        f0(4)*f0(7) + f0(5)*f0(8) + f0(6)*f0(9)
      ee0(6) =        f0(7)*f0(1) + f0(8)*f0(2) + f0(9)*f0(3)

      al0 = 1.d0 - al1

      do i = 1,6
        ee(i) = al0*ee0(i) + al1*ee1(i)
      end do ! i

      end

      subroutine ener3d4(ndm,ndf,nel, d,ul,xl)

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

      integer  a, i,j,l, lint, ndm,ndf,nel
      real*8   detf, xsj, dvol, con1,con2, kineng, poteng
      real*8   fa(3,3),shp(4,8,8), ss(6), dm(6,6)
      real*8   ee(6), sg(4,8), vel(3), mom(3,8), xsm(3), f0(3,3)
      real*8   d(*), ul(ndf,nen,*), xl(ndm,*)

      save

c     Initialize

      do a = 1,nel
        do j = 1,3
          mom(j,a)  = 0.0d0
        end do ! j
      end do ! a

      kineng = 0.0d0
      poteng = 0.0d0

c     Get quadrature information

      if(nint(d(182)).gt.0) then
        call int3dn(nel,lint,sg)
      else
        l = 2
        call int3d(l,lint, sg)
      endif

c     Loop over quadrature points

      do l = 1,lint

        call shp3d(sg(1,l),xsj,shp(1,1,l),xl,ndm,nel)

        dvol = xsj*sg(4,l)

c       Compute deformation gradient and increment

        call fmat3d4(shp(1,1,l),ul,ndf,nel,nen, fa,f0,detf)

c       Compute Green-Lagrange strain

        call emat3d4(1.d0, fa,f0, ee)

        call modl3d4(d,detf,fa, dm,ss,ee,.false.)

c       Accumlate Potential energy for element

        poteng = poteng + (ss(1)*ee(1) + ss(2)*ee(2) + ss(3)*ee(3)
     &                   + ss(4)*ee(4) + ss(5)*ee(5) + ss(6)*ee(6))*dvol

c       Form velocity at gauss point

        con2 = d(4) * dvol
        con1 = con2 * d(7)       ! Consistent factor

c       Inertial contribution from consistent mass approximation

        do i = 1,3

          vel(i) = 0.0d0
          do a = 1,nel
            vel(i) = vel(i) + shp(4,a,l)*ul(i,a,4)
          end do ! a

          do a = 1,nel
            mom(i,a) = mom(i,a) + shp(4,a,l)*vel(i)*con1
          end do ! a

c         Accumulate kinetic energy for element

          kineng = kineng + vel(i)*vel(i)*con1

        end do ! i

c       Inertial contribution from lumped mass approximation

        con2 = con2 - con1       ! Lumped     factor
        do a = 1,nel
          con1 = con2*shp(4,a,l) ! Lumped times shape function
          do i = 1,3
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
        epl(3) = epl(3) + mom(3,a)

c       Angular momentum

        do i = 1,3
          xsm(i) = xl(i,a) + ul(i,a,1)
        end do ! i

        epl(4) = epl(4) + xsm(2)*mom(3,a) - xsm(3)*mom(2,a)
        epl(5) = epl(5) + xsm(3)*mom(1,a) - xsm(1)*mom(3,a)
        epl(6) = epl(6) + xsm(1)*mom(2,a) - xsm(2)*mom(1,a)

      end do ! a

c     Accumulate total kinetic and stored energy

      epl(7) = epl(7) + 0.5d0*kineng
      epl(8) = epl(8) + 0.5d0*poteng

      end

      subroutine fmat3d4(shp,ul,ndf,nel,nen, f,df,detf)

c-----[--.----+----.----+----.-----------------------------------------]
c     Form Deformation Gradient

c     Inputs:

c       shp(4,8)      - shape function derivatives: N_I,i = shp(i,I)
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
      real*8   shp(4,8),ul(ndf,nen,*), f(3,3),df(3,3), detf

      do i = 1,3
        do j = 1,3
          f(i,j)  = 0.0d0
          df(i,j) = 0.0d0
          do k = 1,nel
            f(i,j)  = f(i,j)  + ul(i,k,1)*shp(j,k)
            df(i,j) = df(i,j) + ul(i,k,2)*shp(j,k)
          end do ! k
        end do ! j
        f(i,i) = f(i,i) + 1.0d0
      end do ! i

      detf = f(1,1)*(f(2,2)*f(3,3) - f(2,3)*f(3,2))
     &     + f(1,2)*(f(2,3)*f(3,1) - f(2,1)*f(3,3))
     &     + f(1,3)*(f(2,1)*f(3,2) - f(2,2)*f(3,1))

      end

      subroutine gstf3d4(ss,shp,rhocdv,rholdv,ctan,ndf,nel,nst, s)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Geometric stiffness and consistent mass contribution

c     Input:
c       ss(6)      - Stress times d_volume
c       shp(4,8)   - Shape function derivatives
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

      include 'part0.h'
      include 'rdata.h'

      integer  ndf,nel,nst, a,b, ia, ib, i
      real*8   kgeom, rhocdv,rholdv, lmas,cmas, ctan3
      real*8   ss(6), shp(4,8), s(nst,nst), t(3), ctan(3)

      if(ndfo(1).gt.0 .or. shflg) then
        ctan3 = ctan(3)
      else
        ctan3 = 0.0d0
      endif

      ib = 0
      do b = 1,nel

        t(1) = (ss(1)*shp(1,b)+ss(4)*shp(2,b)+ss(6)*shp(3,b))*ctan(1)
        t(2) = (ss(4)*shp(1,b)+ss(2)*shp(2,b)+ss(5)*shp(3,b))*ctan(1)
        t(3) = (ss(6)*shp(1,b)+ss(5)*shp(2,b)+ss(3)*shp(3,b))*ctan(1)
        lmas = rholdv*shp(4,b)*ctan3
        cmas = rhocdv*shp(4,b)*ctan3

        do i = 1,3
          s(ib+i,ib+i) = s(ib+i,ib+i) + lmas
        end do ! i

        ia = 0
        do a = 1,nel

          kgeom = shp(1,a)*t(1) + shp(2,a)*t(2) + shp(3,a)*t(3)
     &          + shp(4,a)*cmas

          do i = 1,3
            s(ia+i,ib+i) = s(ia+i,ib+i) + kgeom
          end do ! i

          ia = ia + ndf
        end do ! a

        ib = ib + ndf
      end do ! b

      end

      subroutine modl3d4(d,detf,f, dm,ss,ee,dmat)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Finite Deformation Constitutive Models

c     Input:
c         d(*)     - Material parameters
c         f(3,3)   - Deformation gradient
c         detf     - Jacobian determinant at t_n+1

c     Output:
c         dm(6,6)  - Reference (material) elastic moduli
c         ss(6)    - 2nd Piola-Kirchhoff stress tensor
c         ee(6)      Deformation measure
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      logical  dmat
      integer  imat
      real*8   d(*),ss(6),dm(6,6),f(3,3),ee(6), detf

      save

      imat = nint(d(20))
      if(imat.eq.1) then

        call stnh3d4(d,detf,f, dm,ss,dmat)

      elseif(imat.eq.5 .or. imat.eq.6) then

        call stvk3d4(d, ee, dm, ss)

      endif

      end

      subroutine stif3d4(ndm,ndf,nst, d,ul,xl, p,s, isw)

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
      include 'pmod2d.h'
      include 'ptdat6.h'
      include 'tdata.h'

      integer  a,b, i,j,k,l, ia,jb, isw, lint, ndm,ndf,nst

      real*8   fa(3,3),shp(4,8,8), bb(6,3,8), bd(3,6), ss(6), ds(6,6)
      real*8   dvol(8), ee(6), sg(4,8), ace(3), detf, xsj, con1,con2
      real*8   df(3,3), f0(3,3), f1(3,3), bba(6,3,8), bf(3)
      real*8   d(*), ul(ndf,nen,*), xl(ndm,*), p(ndf,*), s(nst,nst)

      save

c     Set body forces

      call sbodyf(d, bf)

c     Get quadrature information

      if(nint(d(182)).gt.0) then
        call int3dn(nel,lint,sg)
      else
        l = 2
        call int3d(l,lint, sg)
      endif

c     Loop over quadrature points

      do l = 1,lint

        call shp3d(sg(1,l),xsj,shp(1,1,l),xl,ndm,nel)

        dvol(l) = xsj*sg(4,l)

c       Compute deformation gradient and increment

        call fmat3d4(shp(1,1,l),ul,ndf,nel,nen, fa,df,detf)

        con1 = (1.d0 - ctan(1))/ctan(1)
        do j = 1,3
          do i = 1,3
            f1(i,j) = fa(i,j) + con1*df(i,j)
            f0(i,j) = fa(i,j) -      df(i,j)
          end do ! i
        end do ! j

c       Compute strain-displacement matrices

        call bmat3d4(fa,shp(1,1,l),bba,nel)
        call bmat3d4(f1,shp(1,1,l),bb ,nel)

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

c       Form acceleration at gauss point

        do i = 1,3
          ace(i) = 0.0d0
          do a = 1,nel
            ace(i) = ace(i) + shp(4,a,l)*ul(i,a,5)
          end do ! a
        end do ! i

c       Compute residual

        con2 = d(4) * dvol(l)
        con1 = con2 * d(7)     ! Consistent factor
        con2 = con2 - con1     ! Lumped     factor
        do a = 1,nel
          do i = 1,3
            do j = 1,6
              p(i,a) = p(i,a) + shp(4,a,l)*dvol(l)*bf(i)
     &                        - bba(j,i,a)*ss(j)
            end do ! j
            p(i,a) = p(i,a) - shp(4,a,l)*(ace(i)*con1 + ul(i,a,5)*con2)
          end do ! i
        end do ! a

c       Compute tangent

        if(isw.eq.3) then

c         Form geometric and inertial parts of tangent stiffness

          if(gflag) then
            call gstf3d4(ss,shp(1,1,l),con1,con2,ctan,ndf,nel,nst, s)
          endif

c         Form material part of tangent

          ia = 0
          do a = 1,nel

            do i = 1,3
              do j = 1,6
                bd(i,j) = 0.0d0
                do k = 1,6
                  bd(i,j) = bd(i,j) + bba(k,i,a)*ds(k,j)
                end do ! k
              end do ! j
            end do ! i

            jb = 0
            do b = 1,nel

              do i = 1,3
                do j = 1,3
                  do k = 1,6
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

      end

      subroutine stnh3d4(d, detf, f, dm, ss, dmat)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: St. Venant - Kirchhoff Material Model
c              Orthotropic elastic material

c     Input:
c       d(*)    - Moduli and Poisson ratios
c       detf    - Determinant of deformation gradient
c       f(9)    - Deformation gradient

c     Outputs:
c       dm(6,6) - Material moduli
c       ss(6)   - 2nd Piola-Kirchhoff stress
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      logical  dmat
      integer  i,j
      real*8   detf, c1,c2, mu, lamb, con1
      real*8   d(*), f(9), dm(6,6), ss(6), cc(6), ci(6)

      save

c     Load moduli

      lamb = d(21)
      mu   = d(22)

c     Compute Green deformation tensor

      cc(1) = f(1)*f(1) + f(2)*f(2) + f(3)*f(3)
      cc(2) = f(4)*f(4) + f(5)*f(5) + f(6)*f(6)
      cc(3) = f(7)*f(7) + f(8)*f(8) + f(9)*f(9)
      cc(4) = f(1)*f(4) + f(2)*f(5) + f(3)*f(6)
      cc(5) = f(4)*f(7) + f(5)*f(8) + f(6)*f(9)
      cc(6) = f(7)*f(1) + f(8)*f(2) + f(9)*f(3)

c     Inverse deformation tensor

      call symi3d4(cc,ci)

c     Compute 2nd P-K stress

      c1  = mu - lamb*detf*(detf - 1.d0)

      do i = 1,6
        ss(i) = - c1*ci(i)
      end do ! i

      do i = 1,3
        ss(i) = ss(i) + mu
      end do ! i

c     Compute tangent tensor

      if(dmat) then

        c1 = 2.d0*c1
        c2 = lamb*detf*(2.d0*detf - 1.d0)

        dm(1,1) = c1*ci(1)*ci(1)
        dm(1,2) = c1*ci(4)*ci(4)
        dm(1,3) = c1*ci(6)*ci(6)
        dm(1,4) = c1*ci(1)*ci(4)
        dm(1,5) = c1*ci(4)*ci(6)
        dm(1,6) = c1*ci(6)*ci(1)

        dm(2,2) = c1*ci(2)*ci(2)
        dm(2,3) = c1*ci(5)*ci(5)
        dm(2,4) = c1*ci(4)*ci(2)
        dm(2,5) = c1*ci(2)*ci(5)
        dm(2,6) = c1*ci(5)*ci(4)

        dm(3,3) = c1*ci(3)*ci(3)
        dm(3,4) = c1*ci(6)*ci(5)
        dm(3,5) = c1*ci(5)*ci(3)
        dm(3,6) = c1*ci(3)*ci(6)

        dm(4,4) = c1*(ci(1)*ci(2) + ci(4)*ci(4))
        dm(4,5) = c1*(ci(4)*ci(5) + ci(2)*ci(6))
        dm(4,6) = c1*(ci(6)*ci(4) + ci(5)*ci(1))

        dm(5,5) = c1*(ci(2)*ci(3) + ci(5)*ci(5))
        dm(5,6) = c1*(ci(5)*ci(6) + ci(3)*ci(4))

        dm(6,6) = c1*(ci(3)*ci(1) + ci(6)*ci(6))

        do i = 1,6
          con1 = c2*ci(i)
          do j = i,6
            dm(i,j) = dm(i,j) + con1*ci(j)
          end do ! j
        end do ! i

c       Multiply shear moduli by half

        do i = 4,6
          do j = i,6
            dm(i,j) = 0.5d0*dm(i,j)
          end do ! j
        end do ! i

c       Make tangent moduli symmetric

        do i = 2,6
          do j = 1,i-1
            dm(i,j) = dm(j,i)
          end do ! j
        end do ! i

      endif

      end

      subroutine stre3d4(ndm,ndf, d,ul,xl, p,s, isw)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Stress output in reference coordinates

c     Input:
c       ndm           - Number of coordinates in mesh per node
c       ndf           - Number of degree-of-freedoms per node
c       nen           - Dimension of stress array
c       d(*)          - Material parameter array
c       ul(ndf,nen,*) - Nodal solution array
c       xl(ndm,*)     - Nodal coordinate array
c       isw           - Switch to control computations

c     Output:
c       p(*)          - Projection weight vector
c       s(nen,*)      - Projection array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'bdata.h'
      include 'cdata.h'
      include 'eldata.h'
      include 'iofile.h'
      include 'strnum.h'

      integer  i,j,k,l, isw, lint, ndm,ndf

      real*8   d(*), p(*), s(nen,*), ul(ndf,*), xl(ndm,*)
      real*8   shp(4,8),shpv(8,8),sg(4,8),df(3,3),detf,xsj,xx1,xx2,xx3
      real*8   f(3,3), ee0(6), ss0(6), ss(6,8), dd(6,6), ee(6,8), spr(3)

      save

c     Get quadrature information and initialize

      if(nint(d(182)).gt.0) then
        call int3dn(nel,lint,sg)
      else
        l = 2
        call int3d(l,lint, sg)
      endif

      do i = 1,6
        ss0(i) = 0.0d0
        ee0(i) = 0.0d0
      end do ! i

      xx1 = 0.0d0
      xx2 = 0.0d0
      xx3 = 0.0d0

c     Loop over quadrature points

      do l = 1,lint

        call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)

c       Compute deformation gradient and increment

        call fmat3d4(shp,ul,ndf,nel,nen, f,df,detf)

c       Compute Green-Lagrange strain

        call emat3d4(1.d0, f,df, ee(1,l))

        call modl3d4(d,detf,f, dd,ss(1,l),ee(1,l),.false.)

        if(isw.eq.4) then

          do i = 1,6
            ss0(i) = ss0(i) + 0.125d0*ss(i,l)
            ee0(i) = ee0(i) + 0.125d0*ee(i,l)
          end do ! i

c         Compute coordinates for outputs

          do i=1,nel
            xsj = 0.125d0*shp(4,i)
            xx1 = xx1 + xsj*xl(1,i)
            xx2 = xx2 + xsj*xl(2,i)
            xx3 = xx3 + xsj*xl(3,i)
          end do ! i

        elseif(isw.eq.8) then

c         Save shape functions for projections

          xsj = xsj*sg(4,l)

          do i = 1,8
            shpv(i,l) = shp(4,i)*xsj
          end do ! i

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

        call pstr3d(ss0,spr)

        write(iow,2002) n,ma,(ss0(j),j=1,6),xx1,xx2,xx3,spr,ee0
        if(ior.lt.0) then
          write(*,2002) n,ma,(ss0(j),j=1,6),xx1,xx2,xx3,spr,ee0
        end if

c     Project stresses

      elseif(isw.eq.8) then

c       Perform integration

        do l = 1,lint

          do i = 1,nel
            p(i) = p(i) + shpv(i,l)
            do k = 1,6
              s(i,k  ) = s(i,k  ) + ss(k,l)*shpv(i,l)
              s(i,k+6) = s(i,k+6) + ee(k,l)*shpv(i,l)
            end do ! k
          end do ! i

        end do ! l

        iste = 12

      end if

2001  format(a1,20a4//5x,'Element Stresses'//
     & '  Elmt  Matl  11-stress  22-stress  33-stress  12-stress',
     &             '  23-stress  31-stress'/
     &             16x,'1-coord    2-coord    3-coord   1-stress',
     &             '   2-stress   3-stress'/
     &           14x,'11-strain  22-strain  33-strain  12-strain',
     &             '  23-strain  31-strain')

2002  format(2i6,1p,6e11.3/12x,1p,6e11.3/12x,1p,6e11.3)

      end

      subroutine stvk3d4(d, ee, dm, ss)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: St. Venant - Kirchhoff Material Model
c              Orthotropic elastic material

c     Input:
c       d(*)    - Moduli and Poisson ratios
c       ee(6)   - Green-Lagrange strains

c     Outputs:
c       dm(6,6) - material moduli
c       ss(6)   - 2nd Piola-Kirchhoff stress
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i,j
      real*8   d(*), dm(6,6), ss(6), ee(6)

c     Load moduli

      do i = 1,6
        do j = 1,6
          dm(i,j) = 0.0d0
        end do ! j
      end do ! i

      do i = 1,3
        dm(i  ,i  ) = d(i+20)
        dm(i+3,i+3) = d(i+26)
        j       = mod(i,3) + 1
        dm(i,j) = d(i+23)
        dm(j,i) = d(i+23)
      end do ! i

c     Compute 2nd P-K stress

      do i = 1,6
        ss(i) = 0.0d0
        do j = 1,6
          ss(i) = ss(i) + dm(i,j)*ee(j)
        end do ! j
      end do ! i

      end

      subroutine symi3d4(cc, ci)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Invert symmetric 3x3 matrix

c     Input:
c       cc(6) - symmetric 3x3 matrix

c     Output:
c       ci(6) - inverse of symmetric 3x3 matrix
c-----[--.----+----.----+----.-----------------------------------------]

      implicit none

      integer  i
      real*8   cc(6), ci(6), rdetc

      ci(1) = cc(2)*cc(3) - cc(5)*cc(5)
      ci(2) = cc(3)*cc(1) - cc(6)*cc(6)
      ci(3) = cc(1)*cc(2) - cc(4)*cc(4)
      ci(4) = cc(5)*cc(6) - cc(3)*cc(4)
      ci(5) = cc(6)*cc(4) - cc(1)*cc(5)
      ci(6) = cc(4)*cc(5) - cc(2)*cc(6)

      rdetc = 1.d0/(cc(1)*ci(1) + cc(4)*ci(4) + cc(6)*ci(6))
      do i = 1,6
        ci(i) = ci(i)*rdetc
      end do ! i

      end
