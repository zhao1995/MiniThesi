c$Id:$
      subroutine mnrv3f(d,detf,bb,gg, sig,aa,xlamd,ha,estore)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension detf(3), use in fengy                  22/11/2011
c       2. Correct computation of stress and tangent from   30/04/2013
c          second invariant term
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Finite deformation elasticity Mooney-Rivlin model
c         W(I_1,I_2,J) = U(J) + 0.5*mu*{(1-c)*[I_1 - 3 - 2 log(J)]
c                                        + c *[I_2 - 3 - 4 log(J)]
c                   I_1 = bb:1  (trace b)
c                   I_2 = 0.5*(I_1*I_1 - b:b)

c         Note: Coding in terms of F - I = G

c     Inputs:
c         d(100)     Material parameters
c         d(21)      Bulk  modulus
c         d(22)      Shear modulus
c         detf(3)    Jacobian determinant at t_n+1
c         bb(6)      Left Cauchy-Green deformation tensor
c         gg(6)      Left Cauchy-Green displacement gradient tensor

c     Outputs:
c         sig(6)     CAUCHY stress tensor
c         aa(6,6)    CAUCHY (spatial) elastic moduli
c         estore     Stored energy density
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,k,l, ii,jj, jsw
      real*8    detf(3), press,u,up,upp, mu, mu2, ha,hp,hpp, estore
      real*8    xlamd,logj, nu,nu2,nu4, trg, trg2
      real*8    d(*),sig(6),aa(6,6),bb(6),gg(6), g2(6)
      real*8    b(3,3),c(9,9)
      integer   ia(9),ja(9)

      data      ia / 1, 2, 3, 1, 2, 2, 3, 3, 1 /
      data      ja / 1, 2, 3, 2, 1, 3, 2, 1, 3 /

c     Compute pressure and its derivative

      jsw = nint(d(170))
      call fengy3(d,detf(3),u,up,upp, ha,hp,hpp, jsw)
      press =  up  + xlamd * hp
      upp   = (upp + xlamd * hpp) * detf(1)

c     Set CAUCHY stresses and elastic tangent moduli

      nu  =  d(22)*d(23)/detf(1)
      mu  =  d(22)/detf(1) - nu
      nu2 =  nu + nu
      nu4 =  nu2 + nu2
      mu2 =  mu + mu

c     First invariant term

      do i = 1,3
        sig(i  ) = mu*gg(i) + press
        sig(i+3) = mu*gg(i+3)
      end do ! i

c     Square of b tensor [omits 1.]

      g2(1) = gg(1)*gg(1) + gg(4)*gg(4) + gg(6)*gg(6) + 2.d0*gg(1)
      g2(2) = gg(4)*gg(4) + gg(2)*gg(2) + gg(5)*gg(5) + 2.d0*gg(2)
      g2(3) = gg(6)*gg(6) + gg(5)*gg(5) + gg(3)*gg(3) + 2.d0*gg(3)
      g2(4) = gg(1)*gg(4) + gg(4)*gg(2) + gg(6)*gg(5) + 2.d0*gg(4)
      g2(5) = gg(4)*gg(6) + gg(2)*gg(5) + gg(5)*gg(3) + 2.d0*gg(5)
      g2(6) = gg(6)*gg(1) + gg(5)*gg(4) + gg(3)*gg(6) + 2.d0*gg(6)

c     Trace term [omits 3.0 term]

      trg   = gg(1) + gg(2) + gg(3)

c     Set stress

      do i = 1,3
        sig(i  ) = sig(i  ) + nu*(3.d0*gg(i) + trg + trg*gg(i) - g2(i))
        sig(i+3) = sig(i+3) + nu*(3.d0*gg(i+3) + trg*gg(i+3) - g2(i+3))
      end do ! i

c     Tangent terms from Kronnecker deltas of I_1 and I_2

      do i = 1,3
        aa(i  ,i  ) = mu2 + nu4 - press + upp
        aa(i+3,i+3) = mu  + nu2 - press
      end do ! i

c     Tangent terms from second invariant

      b(1,1) = bb(1)
      b(1,2) = bb(4)
      b(1,3) = bb(6)
      b(2,1) = bb(4)
      b(2,2) = bb(2)
      b(2,3) = bb(5)
      b(3,1) = bb(6)
      b(3,2) = bb(5)
      b(3,3) = bb(3)

      c = 0.d0
      do jj = 1,9
        k = ia(jj)
        l = ja(jj)
        do ii = 1,9
          i = ia(ii)
          j = ja(ii)
          c(ii,jj) = c(ii,jj) + b(i,j)*b(k,l)
     &                 - 0.5d0*(b(i,k)*b(j,l) + b(i,l)*b(j,k))
        end do ! ii
      end do ! jj

      c(:,4) = 0.5d0*(c(:,4) + c(:,5))
      c(:,5) = 0.5d0*(c(:,6) + c(:,7))
      c(:,6) = 0.5d0*(c(:,9) + c(:,9))
      c(4,:) = 0.5d0*(c(4,:) + c(5,:))
      c(5,:) = 0.5d0*(c(6,:) + c(7,:))
      c(6,:) = 0.5d0*(c(8,:) + c(9,:))

      aa(1:6,1:6) = aa(1:6,1:6) + nu2*c(1:6,1:6)

c     Add volumetric correction to aa

      upp     = press   + upp
      aa(1,2) = aa(1,2) + upp
      aa(2,1) = aa(1,2)
      aa(1,3) = aa(1,3) + upp
      aa(3,1) = aa(1,3)
      aa(2,3) = aa(2,3) + upp
      aa(3,2) = aa(2,3)

c     Compute stored energy density

      logj   =  log(abs(detf(1)))
      trg2   =  g2(1)*g2(1) + g2(2)*g2(2) + g2(3)*g2(3)
     &       + (g2(4)*g2(4) + g2(5)*g2(5) + g2(6)*g2(6))*2.d0
      estore =  u + 0.5d0*d(22)*((1.d0 - d(23))*(trg - 2.0d0*logj)
     &       + d(23)*(2.d0*trg + 0.5d0*(trg*trg - trg2) - 4.0d0*logj))

      end
