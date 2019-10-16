c$Id:$
      function ylds2d(d,sigtr,alp,eta,dd,beta,epn,lam0,hprim)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c     2. Replace sqrt2 by stwo                              18/03/2010
c-----[--.----+----.----+----.-----------------------------------------]
c     Plane stress plasticity routine for return map algorithm
c     Inputs:
c        d(*)     - Parameters
c        sigtr(*) - Trial stress
c        alp(*)   - Back stress
c        eta(*)   - Plastic strain

c     Outputs:
c        dd(6,6)  - Moduli
c        beta     - Yield measure
c        ylds2d   - Yield function value
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pconstant.h'
      include  'rdata.h'

      integer   i,icnt
      real*8    d(*),sigtr(3),alp(3),eta(3),psi(3),dd(6,6)
      real*8    toli,stwo,psi1,psi2,twoh3,d1,d2,e2,e3
      real*8    f1,f2,ff1,ff2,phit,ep,rad2,f3,phi,dphi,dlam,lam
      real*8    gam1,gam2,beta,lam0,epn,hprim,xkprim, ylds2d,hard2d

      save

      data      toli  /1.d-8/

      stwo   = sqrt2*0.5d0
      psi(1) = stwo*( eta(1) + eta(2))
      psi(2) = stwo*(-eta(1) + eta(2))
      psi(3) = eta(3)

c     Compute scale factor

      psi1 = psi(1)*psi(1)*one3
      psi2 = psi(2)*psi(2) + psi(3)*psi(3)*2.d0
      twoh3= two3*hprim
      d1   = d(1)/(1.d0 - d(2))*one3
      d2   = d(1)/(1.d0 + d(2))
      e3   = d1 + twoh3
      e2   = d2 + twoh3

c     Newton's method for determining correct lambda

      lam  = lam0
      icnt = 0
100   f1   = 1.d0/(1.d0 + e3*lam)
      f2   = 1.d0/(1.d0 + e2*lam)
      ff1  = f1*f1*psi1
      ff2  = f2*f2*psi2
      phit = 0.50d0*(ff1  + ff2)
      ep   = epn + sqrt(four3*phit)*lam
      rad2 = hard2d(d,ep,xkprim,hprim)
      f3   = (1.d0 - lam*two3*xkprim)
      phi  = phit - rad2
      dphi = f3*(e3*ff1*f1 + e2*ff2*f2) + four3*xkprim*phit
      dlam = phi/dphi
      if(dlam.gt.0.0d0.or.abs(dlam).lt.abs(lam)) then
        lam = lam + dlam
      else
        lam = 0.50d0*abs(lam)
      endif
      icnt = icnt + 1
      if(icnt.gt.100) go to 110
      if(abs(dlam).gt.toli*abs(lam) .and. rnmax.gt.0.0d0) go to 100
      go to 120
110   write(iow,2000) dlam,lam

c     Scale psi onto yield surface using lambda

120   continue

c     Compute plasticity map dd array

      f1      = 1.d0/(1.d0 + e3*lam)
      f2      = 1.d0/(1.d0 + e2*lam)
      gam1    = 1.d0 + twoh3*lam
      gam2    = 1.d0 - two3*xkprim*lam
      d1      = 1.5d0*d1*f1*gam1
      d2      = 0.5d0*d2*f2*gam1

      dd(1,1) = d1 + d2
      dd(1,2) = d1 - d2
      dd(2,1) = d1 - d2
      dd(2,2) = d1 + d2
      dd(4,4) = d2

c     Compute stresses on yield surface

      psi(1)  = psi(1)*f1
      psi(2)  = psi(2)*f2
      psi(3)  = psi(3)*f2
      eta(1)  = stwo*(psi(1) - psi(2))
      eta(2)  = stwo*(psi(1) + psi(2))
      eta(3)  = psi(3)
      ylds2d  = (eta(1)**2 + eta(2)**2 - eta(1)*eta(2))*one3
     &        +  eta(3)**2
      rad2    = hard2d(d,ep,xkprim,hprim)
      do i = 1,3
        alp(i)  = alp(i) + lam*twoh3*eta(i)
        sigtr(i)= eta(i) + alp(i)
      end do ! i
      beta   = four3*ylds2d*(gam1*xkprim + gam2*hprim)*gam1/gam2
      ylds2d = sqrt(ylds2d*3.d0)/d(41)
      lam0   = lam
      epn    = ep

2000  format(' * * Warning * * Failure to converge in ylds2d'/
     &       '     dlam =',1p,1e12.5,' lam =',1p,1e12.5)

      end
