c$Id:$
      subroutine stvktm(d, detf, f, gradt,tg,
     &                  sig,flux, dd,kt, energy)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add .false. to call of tranr4                    06/12/2011
c       2. Correct initial definitions for ee(1:3) (2*)     02/04/2012
c       3. Push stress with full deformation gradient       16/07/2013
c          Add by pushf2 instead of pushr2.
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: St. Venant-Kirchhoff Material: Orthotropic

c     Input:
c       d(*)     - Moduli and Poisson ratios
c       f(9,2)   - Deformation gradient at time: t_n+1, t_n
c       gradt(3) - Spatial thermal gradient
c       tg       - Temperature

c     Outputs:
c       sig(1:6) - Cauchy stress                           : (momentum)
c       flux(1:3)- Thermal flux                            : (heat)
c       flux(4)  - Thermal expansion and structural heating: (heat)
c                 (on current configuration)

c       dd(7,7)  - Spatial moduli                          : (momentum)
c                  1:6,7 = Thermal    coupling             : (momentum)
c                  7,1:6 = Mechanical coupling             : (heat)
c       kt(3,3)  - Thermal conductivity                    : (momentum)

c       energy   - Energy density
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'ddata.h'
      include   'tdata.h'

      real*8     detf, tg, energy, R
      real*8     d(*),  f(9,2),  gradt(3)
      real*8     sig(6),flux(4), dd(7,7), kt(3,3)

      integer    i,j
      real*8     delt, tgdt
      real*8     dm(6,6), ss(6) , ee(6), de(6), alp(3), b0(6),bhat(6)
      real*8     ds(6,6), t(6,6), k0(3,3)

      save

c     MOMENTUM EQUATION TERMS

c     Set reference configuration moduli

      do j = 1,6
        do i = 1,6
          dm(i,j) = 0.0d0
        end do ! i
      end do ! j

      do i = 1,3

c       Diagonal terms

        dm(i  ,i  ) = d(i+20)
        dm(i+3,i+3) = d(i+26)

c       Off-diagonal terms

        j       = mod(i,3) + 1
        dm(i,j) = d(i+23)
        dm(j,i) = d(i+23)
      end do ! i

c     Compute thermal expansion terms

      do i = 1,3
        alp(i)   = d(46+i)
      end do ! i

      do i = 1,6
        b0(i) = dm(i,1)*alp(1) + dm(i,2)*alp(2) + dm(i,3)*alp(3)
      end do ! i

c     Compute Green-Lagrange strain

      ee(1) = f(1,1)*f(1,1) + f(2,1)*f(2,1) + f(3,1)*f(3,1)
     &      + f(1,1)*2.d0
      ee(2) = f(4,1)*f(4,1) + f(5,1)*f(5,1) + f(6,1)*f(6,1)
     &      + f(5,1)*2.d0
      ee(3) = f(7,1)*f(7,1) + f(8,1)*f(8,1) + f(9,1)*f(9,1)
     &      + f(9,1)*2.d0

      ee(4) = f(1,1)*f(4,1) + f(2,1)*f(5,1) + f(3,1)*f(6,1)
     &      + f(4,1) + f(2,1)
      ee(5) = f(4,1)*f(7,1) + f(5,1)*f(8,1) + f(6,1)*f(9,1)
     &      + f(8,1) + f(6,1)
      ee(6) = f(7,1)*f(1,1) + f(8,1)*f(2,1) + f(9,1)*f(3,1)
     &      + f(3,1) + f(7,1)

c     Compute increment in Green-Lagrange strain: dE = (E_n+1 - E_n)

      de(1) = ee(1) - (f(1,2)*f(1,2) + f(2,2)*f(2,2) + f(3,2)*f(3,2))
     &      - f(1,2)*2.d0
      de(2) = ee(2) - (f(4,2)*f(4,2) + f(5,2)*f(5,2) + f(6,2)*f(6,2))
     &      - f(5,2)*2.d0
      de(3) = ee(3) - (f(7,2)*f(7,2) + f(8,2)*f(8,2) + f(9,2)*f(9,2))
     &      - f(9,2)*2.d0

      de(4) = ee(4) - (f(1,2)*f(4,2) + f(2,2)*f(5,2) + f(3,2)*f(6,2))
     &      - f(4,2) + f(2,2)
      de(5) = ee(5) - (f(4,2)*f(7,2) + f(5,2)*f(8,2) + f(6,2)*f(9,2))
     &      - f(8,2) + f(6,2)
      de(6) = ee(6) - (f(7,2)*f(1,2) + f(8,2)*f(2,2) + f(9,2)*f(3,2))
     &      - f(3,2) + f(7,2)

c     Modify Green-Lagrange strain for thermal strain and half

      delt = tg - d(9)
      do i = 1,3
        ee(i) = ee(i)*0.5d0 - alp(i)*delt
        de(i) = de(i)*0.5d0
      end do ! i

c     Compute 2nd P-K stress

      do i = 1,6
        ss(i) = 0.0d0
        do j = 1,6
          ss(i) = ss(i) + dm(i,j)*ee(j)
        end do ! j
      end do ! i

c     Push to current configuration: Store thermal in col 7 of dd array

      call tranr4(f,f,t,.true.)
      call pushr4(t,t,dm, ds,detf)
      call pushf2(f,ss,sig,detf)

      do j = 1,6
        do i = 1,6
          dd(i,j) = ds(i,j)
        end do ! i
      end do ! j

c     Thermal stress compling tangent

      call pushr2(f,b0,bhat,detf)
      do i = 1,6
        dd(i,7) = -bhat(i)
      end do ! i

c     THERMAL PROPERTIES

c     Set reference conductivity

      do j = 1,3
        do i = 1,3
          k0(i,j) = 0.0d0
        end do ! i
        k0(j,j) = d(60+j)
      end do ! j

c     Push to current configuration (divide by det F)

      call pusht2(f,k0,kt,detf)

c     Form flux in current configuration

      do i = 1,3
        flux(i) = -(kt(i,1)*gradt(1)
     &            + kt(i,2)*gradt(2)
     &            + kt(i,3)*gradt(3))
      end do ! i

c     Compute thermal expansion heat: (dot(E) * D * alpha)/J
c     (Push heating term to current configuration)

      if(dt.gt.0.0d0) then
        tgdt = (tg + d(129))/dt
        do i = 1,6
          dd(7,i) = bhat(i)*tgdt/detf
        end do ! i

        R = (de(1)*b0(1) + de(2)*b0(2) + de(3)*b0(3)
     &    +  de(4)*b0(4) + de(5)*b0(5) + de(6)*b0(6))/dt
      else
        do i = 1,6
          dd(7,i) = 0.0d0
        end do ! i
        R = 0.0d0
      endif
      dd(7,7) = R/detf
      flux(4) = dd(7,7)*(tg + d(129))

c     Compute energy density

      energy = 0.5d0*(ss(1)*ee(1) + ss(2)*ee(2) + ss(3)*ee(3)
     &              + ss(4)*ee(4) + ss(5)*ee(5) + ss(6)*ee(6))

      end
