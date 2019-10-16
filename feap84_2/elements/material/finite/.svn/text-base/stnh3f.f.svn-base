c$Id:$
      subroutine stnh3f(d,detf, db, sig,aa,xlamd,ha,estore)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension detf(3), use in fengy                  22/11/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Finite deformation elasticity Neo-Hookean model

c     Inputs:
c         d(100)     Material parameters
c         d(21)      Bulk  modulus
c         d(22)      Shear modulus
c         detf(1)    Jacobian determinant at t_n+1
c         detf(3)    Jacobian determinant at t_n+1 - 1.0d0
c         g(3,3)     Displacement gradient
c         db(6)      Left Cauchy-Green tensor - 1.0d0

c     Outputs:
c         sig(6)     CAUCHY stress tensor
c         aa(6,6)    CAUCHY (spatial) elastic moduli
c         estore     Stored energy density
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'debugs.h'

      integer   i,jsw
      real*8    detf(3), press,u,up,upp, mu, mu2, ha,hp,hpp, estore
      real*8    xlamd,logj
      real*8    d(*),sig(6),aa(6,6),db(6)

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
        aa(i  ,i  ) = mu2 - press + upp
        aa(i+3,i+3) = mu  - press
      end do ! i

c     Add volumetric correction to aa

      upp     = press   + upp
      aa(1,2) = aa(1,2) + upp
      aa(2,1) = aa(1,2)
      aa(1,3) = aa(1,3) + upp
      aa(3,1) = aa(1,3)
      aa(2,3) = aa(2,3) + upp
      aa(3,2) = aa(2,3)

c     Compute stored energy

      logj   = log(abs(detf(1)))
      estore = u + d(22)*(0.5d0*(db(1) + db(2) + db(3)) - logj)

      end
