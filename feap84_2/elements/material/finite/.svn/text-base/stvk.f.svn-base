c$Id:$
      subroutine stvk(d, detf, fa,df, sig,ds, energy)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Always compute energy just before exit           25/04/2011
c       2. Add flag to call of tranr4                       06/12/2011
c       3. Add rotation vector option for moduli            02/05/2013
c       4. Push stress with full deformation gradient       16/07/2013
c          Add by pushf2 instead of pushr2.
c       5. Change outputs from 'print' to 'mprint'          31/10/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: St. Venant-Kirchhoff Material: Orthotropic

c     Input:
c       d(*)    - Moduli and Poisson ratios
c       fa(9)   - Deformation gradient at time: t_n+a
c       df(9)   - Incremental Deformation gradient at time: t_n+a

c     Outputs:
c       sig(6)  - Cauchy stress
c       ds(6,6) - Spatial moduli
c       energy  - Energy density
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'ddata.h'
      include 'debugs.h'

      integer  i,j
      real*8   detf, fac1,fac2, energy, dot
      real*8   d(*), fa(9), df(9),  dm(6,6), ss(6), ee(6)
      real*8   fn(9),f1(9), sig(6), ds(6,6), tl(6,6),tr(6,6)
      real*8   ro(3,3),tt(6,6),dd(6,6)

      save

c     Set moduli

      do j = 1,6
        do i = 1,6
          dd(i,j) = 0.0d0
        end do ! i
      end do ! j
      do i = 1,3
        dd(i  ,i  ) = d(i+20)
        dd(i+3,i+3) = d(i+26)
      end do

      do i = 1,3
        j       = mod(i,3) + 1
        dd(i,j) = d(i+23)
        dd(j,i) = d(i+23)
      end do ! i

c     Rotate moduli for specified vectors

      if(nint(d(242)).eq.1) then
        do i = 1,3
          ro(i,1) = d(242+i)
          ro(i,2) = d(245+i)
        end do ! i
        call triad(ro)
        call tranr4(ro,ro,tt,.false.)
        call pushr4(tt,tt,dd,dm,1.d0)
      else
        dm = dd
      endif

c     Compute factors for updates

      if(energy.ge.0.0d0) then
        fac1 = 0.5d0 * theta(3)
        fac2 = 1.d0/theta(3) - 1.d0
      else
        fac1 = 0.5d0
        fac2 = 0.0d0
      endif

c     Compute deformation gradients at t_n and t_n+1

      do i = 1,9
        fn(i) = fa(i) -      df(i)
        f1(i) = fa(i) + fac2*df(i)
      end do ! i

c     Compute Green-Lagrange strains

      fac2  = 0.5d0 - fac1

      ee(1) = fac1*(f1(1)*f1(1) + f1(2)*f1(2) + f1(3)*f1(3))
     &      + fac2*(fn(1)*fn(1) + fn(2)*fn(2) + fn(3)*fn(3))
     &      + 2.d0*(fac1*f1(1) + fac2*fn(1))
      ee(2) = fac1*(f1(4)*f1(4) + f1(5)*f1(5) + f1(6)*f1(6))
     &      + fac2*(fn(4)*fn(4) + fn(5)*fn(5) + fn(6)*fn(6))
     &      + 2.d0*(fac1*f1(5) + fac2*fn(5))
      ee(3) = fac1*(f1(7)*f1(7) + f1(8)*f1(8) + f1(9)*f1(9))
     &      + fac2*(fn(7)*fn(7) + fn(8)*fn(8) + fn(9)*fn(9))
     &      + 2.d0*(fac1*f1(9) + fac2*fn(9))

      fac1  = fac1 + fac1
      fac2  = fac2 + fac2

      ee(4) = fac1*(f1(1)*f1(4) + f1(2)*f1(5) + f1(3)*f1(6))
     &      + fac2*(fn(1)*fn(4) + fn(2)*fn(5) + fn(3)*fn(6))
     &      + fac1*(f1(4) + f1(2)) + fac2*(fn(4) + fn(2))
      ee(5) = fac1*(f1(4)*f1(7) + f1(5)*f1(8) + f1(6)*f1(9))
     &      + fac2*(fn(4)*fn(7) + fn(5)*fn(8) + fn(6)*fn(9))
     &      + fac1*(f1(8) + f1(6)) + fac2*(fn(8) + fn(6))
      ee(6) = fac1*(f1(7)*f1(1) + f1(8)*f1(2) + f1(9)*f1(3))
     &      + fac2*(fn(7)*fn(1) + fn(8)*fn(2) + fn(9)*fn(3))
     &      + fac1*(f1(7) + f1(3)) + fac2*(fn(7) + fn(3))

c     Compute 2nd P-K stress

      do i = 1,6
        ss(i) = 0.0d0
        do j = 1,6
          ss(i) = ss(i) + dm(i,j)*ee(j)
        end do ! j
      end do ! i

c     Push to current configuration

      if(energy.ge.0.0d0) then

        call tranr4(fa,fa,tl,.true.)
        call tranr4(f1,fa,tr,.true.)
        call pushr4(tl,tr,dm, ds,detf)
        call pushf2(fa,ss,sig,detf) ! adds identity to fa

      endif

c     Compute energy density

      energy = 0.5d0*dot(ss,ee,6)

      if(debug) then
        call mprint( ee,1,6,1,'EE')
        call mprint( ss,1,6,1,'SS')
        call mprint(sig,1,6,1,'SIG')
      endif

      end
