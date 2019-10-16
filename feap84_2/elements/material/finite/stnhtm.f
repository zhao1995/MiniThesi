c$Id:$
      subroutine stnhtm(d,detf,f, db, gradt, tg, sig,ds, flux,kt,
     &                  xlamd,ha,estore)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension detf(3), use in fengy                  22/11/2011
c       2. Pass detf(1) to pushr2                           09/01/2012
c       3. Compute thermal expansion from sum of alphas in  10/07/2012
c          d(3) -- originally d(3) only had alp(1)
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Finite deformation elasticity Neo-Hookean model

c     Inputs:
c         d(100)     Material parameters
c         d(21)      Bulk  modulus
c         d(22)      Shear modulus
c         detf(1)    Jacobian determinant at t_n+1
c         detf(3)    Jacobian determinant at t_n+1 - 1.0d0
c         f(3,3)     Deformation gradient
c         db(6)      Left Cauchy-Green tensor - 1.0d0
c         gradt(3)   Gradient of temperature (spatial)
c         tg         Temperature

c     Outputs:
c         sig(6)     CAUCHY stress tensor
c         flux(4)    Thermal spatial flux
c         ds(7,7)    CAUCHY (spatial) elastic moduli
c         kt(3,3)    Thermal conductivity
c         estore     Stored energy density
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'debugs.h'
      include  'pconstant.h'

      integer   i,j, jsw
      real*8    detf(3), press,u,up,upp, mu, mu2, ha,hp,hpp, estore
      real*8    tg, xlamd,logj, sigtg, delt
      real*8    d(*),sig(6),ds(7,7),f(3,3),db(6), sigt(6), sigt0(6)
      real*8    gradt(3), flux(4), kt(3,3), k0(3,3)

c     Compute pressure and its derivative

      jsw = nint(d(170))
      call fengy3(d,detf(3),u,up,upp, ha,hp,hpp, jsw)
      press =  up  + xlamd * hp
      upp   = (upp + xlamd * hpp) * detf(1)

c     Set CAUCHY stresses and elastic tangent moduli

      mu  =  d(22)/detf(1)
      mu2 =  mu + mu
      do i = 1,3
        sig(i  )    = mu*db(i) + press
        sig(i+3)    = mu*db(i+3)
        ds(i  ,i  ) = mu2 - press + upp
        ds(i+3,i+3) = mu  - press
      end do ! i

c     Add volumetric correction to ds

      upp     = press   + upp
      ds(1,2) = ds(1,2) + upp
      ds(2,1) = ds(1,2)
      ds(1,3) = ds(1,3) + upp
      ds(3,1) = ds(1,3)
      ds(2,3) = ds(2,3) + upp
      ds(3,2) = ds(2,3)

c     Thermal strain modifications

      sigtg = (d(21) + two3*d(22))*d(3)

      do i = 1,3
        sigt0(i)   = sigtg
        sigt0(i+3) = 0.0d0
      end do ! i

      call pushr2(f,sigt0,sigt,detf(1))
      delt  =  tg - d(9)
      do i = 1,6
        sig(i)  =  sig (i) - sigt(i)*delt
        ds(i,7) = -sigt(i)
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

      call pusht2(f,k0,kt,detf(1))

c     Form flux in current configuration

      do i = 1,3
        flux(i) = -(kt(i,1)*gradt(1)
     &            + kt(i,2)*gradt(2)
     &            + kt(i,3)*gradt(3))
      end do ! i

c     Compute stored energy

      logj   = log(abs(detf(1)))
      estore = u + d(22)*(0.5d0*(db(1) + db(2) + db(3)) - logj)

      end
