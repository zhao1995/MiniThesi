      subroutine pengy3d(ndm,ndf,nel, weng, d,ul,xl)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Momentum/energy computation

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

      integer  a, i,l, ndm,ndf,nel
      real*8   con1,con2, kineng, poteng
      real*8   d(*), ul(ndf,nen,*), xl(ndm,*), weng(*)
      real*8   vel(3), mom(3,nel), xsm(3)

      save

c     Initialize

      mom(:,:) = 0.0d0
      kineng   = 0.0d0
      poteng   = 0.0d0

c     Loop over quadrature points

      do l = 1,lint

c       Accumlate Potential energy for element

        poteng = poteng + weng(l)*jac(l)

c       Form velocity at gauss point

        con2 = d(4) * jac(l)
        con1 = con2 * d(7)       ! Consistent factor

c       Inertial contribution from consistent mass approximation

        vel(:) = 0.0d0
        do a = 1,nel
          vel(:) = vel(:) + shp3(4,a,l)*ul(:,a,4)
        end do ! a

        do a = 1,nel
          mom(:,a) = mom(:,a) + shp3(4,a,l)*vel(:)*con1
        end do ! a

c       Accumulate kinetic energy for element

        do i = 1,ndm
          kineng = kineng + vel(i)**2*con1
        end do ! i

c       Inertial contribution from lumped mass approximation

        con2 = con2 - con1       ! Lumped     factor
        do a = 1,nel
          con1 = con2*shp3(4,a,l) ! Lumped times shape function
          mom(:,a) = mom(:,a) + ul(:,a,4)*con1
          do i = 1,3
            kineng   = kineng + ul(i,a,4)**2*con1
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

      epl(7) = epl(7) + kineng*0.5d0
      epl(8) = epl(8) + poteng

      end
