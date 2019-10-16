c$Id:$
      subroutine estrsd(d,ta,eps,sig,dd,dr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Linear Elastic Constitutive Model

c     Inputs:
c        d(*)      - Parameters
c        ta        - Temperature change from stress free state
c        eps(*)    - Strain

c     Outputs:
c        sig(*)    - Stress
c        dd(*,*,2) - Moduli
c        dr(*,*)   - Moduli for Rayleigh damping
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      include  'elcoor.h'
      include  'elengy.h'

      integer   i,j
      real*8    ta,psi, d(*),eps(6),sig(6),dd(6,6,2),dr(6,6),beta(6)
      real*8    tau(6), tt(6,6), dm(6,6),betm(6), ro(3,3)

      save

c     Extract orthotropic vectors

      if(nint(d(242)).eq.1) then
        do i = 1,3
          ro(i,1) = d(242+i)
          ro(i,2) = d(245+i)
        end do ! i
        call triad(ro)
        psi = 0.0d0
        call dmat2d(d,psi,dm,betm)      ! Transform thermal stress
        call pushr2(ro,betm,beta,1.0d0)
        call tranr4(ro,ro,tt,.false.)
        call pushr4(tt,tt,dd,dm,1.d0)   ! Transform moduli

c     Set rotation angle in 1-2 plane

      else
        if(nint(d(85)).gt.0) then
          psi = atan2(xref(2)-d(87),xref(1)-d(86))
        else
          psi = d(31)
        endif
        call dmat2d(d,psi,dd,beta)
      endif

c     Compute stress

      do i = 1,6
        sig(i)    = sig(i) - beta(i)*ta       ! Thermal stress
        tau(i)    = 0.0d0
        dd(i,1,2) = -beta(i)
        do j = 1,6
          tau(i)  = tau(i) + dd(i,j,1)*eps(j) ! Mechanical stress
          dr(i,j) = dd(i,j,1)                 ! Rayleigh damping moduli
        end do ! j
      end do ! i

c     Compute stored energy and final stress

      estore = 0.0d0
      do i = 1,6
        estore = estore + (sig(i) + 0.5d0*tau(i))*eps(i)
        sig(i) = sig(i) + tau(i)
      end do ! i

c     Set plane stress case (dd(3,3) = 0.0d0)

      if(dd(3,3,1) .eq. 0.0d0 ) then
        eps(3) = d(90)*sig(1) + d(91)*sig(2) + d(92)*ta
        sig(3) = 0.0d0
      endif

      end
