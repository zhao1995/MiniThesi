c$Id:$
      subroutine modmnrv(d,f,finv,detf,b1, hn,ntm, sig,aa,xlamd,ha,engy)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'pconstant.h'                 14/11/2006
c       2. Add viscoelastic option and correct definition   07/05/2007
c          of energy in comment below.
c       3. Dimension detf(3), use in fengy                  22/11/2011
c       4. Pass detf(1) to viscfd                           01/05/2012
c       5. Change bd to bb in stress computation            23/05/2013
c          Modify last term of tangent for 9x9 to 6x6 form
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compressible Mooney-Rivlin model
c                (J_2/3 regularization)
c                _ __      _                   _             _
c              W(J,bb) = U(J) + 0.5*mu*[(1-c)*(I_1 - 3) + c*(I_2 - 3)]
c                  __
c                  bb  = J^(-2/3) * b1
c                  _     __                    __
c                  I_1 = bb:1           (trace bb)
c                  _          _   _     __ __
c                  I_2 = 0.5*(I_1*I_1 - bb:bb)

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'
      include  'elcount.h'
      include  'elengy.h'
      include  'iofile.h'

      integer   i,j,k,l, ii,jj, nstv, ntm, jsw, ia(9),ja(9)
      integer   nprop, nvect
      real*8    detf(3), detfi, j23, trbb3, bdi, mub1,mub2,mub3
      real*8    nub1,nub2,nub3, i1b,i2b
      real*8    u, up, upp, ha, hp, hpp, press, xlamd, engy
      real*8    d(*), f(3,3), finv(3,3), b1(6), hn(*), sig(6), aa(6,6)
      real*8    bb(6),bd(6),b2(6), b2d(6), dd(9,9), b(3,3)
      real*8    de(3,3)

      save

      data      ia / 1, 2, 3, 1, 2, 2, 3, 3, 1 /
      data      ja / 1, 2, 3, 2, 1, 3, 2, 1, 3 /
      data      de / 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0 /

c     Compute deviatoric be

      detfi = 1.d0/detf(1)
      j23   = detfi**two3
      do i = 1,6
        bb(i) = b1(i) * j23
      end do ! i

c     Compute bb:bb

      b2(1) = bb(1)*bb(1) + bb(4)*bb(4) + bb(6)*bb(6)
      b2(2) = bb(4)*bb(4) + bb(2)*bb(2) + bb(5)*bb(5)
      b2(3) = bb(6)*bb(6) + bb(5)*bb(5) + bb(3)*bb(3)
      b2(4) = bb(1)*bb(4) + bb(4)*bb(2) + bb(6)*bb(5)
      b2(5) = bb(4)*bb(6) + bb(2)*bb(5) + bb(5)*bb(3)
      b2(6) = bb(6)*bb(1) + bb(5)*bb(4) + bb(3)*bb(6)

      i1b   = bb(1) + bb(2) + bb(3)
      i2b   = (i1b*i1b - b2(1) - b2(2) - b2(3))*0.5d0

c     First invariant strain

      trbb3 = i1b   * one3
      do i = 1,3
        bd(i  ) = bb(i  ) - trbb3
        bd(i+3) = bb(i+3)
      end do ! i

c     Second invariant strain

      do i = 1,3
        b2d(i  ) = i1b*bb(i  ) - b2(i  ) - two3*i2b
        b2d(i+3) = i1b*bb(i+3) - b2(i+3)
      end do ! i

c     Compute deviatoric Kirchhoff stress tensor.

      nub1 = d(22)*d(23)            ! mu*c
      mub1 = d(22) - nub1           ! mu*(1 - c)
      do i = 1,6
        sig(i) = mub1 * bd(i) + nub1*b2d(i)
      end do ! i

c     Compute tangent tensor
c                                  __             __     _
c     Rank one update: -2/3 mu * ( bd x g +  g x  bd ) / J

      mub3 = two3 * mub1
      do i = 1,6
        bdi = mub3 * bd(i)
        do j = 1,3
          aa(i,j) =  aa(i,j) - bdi
          aa(j,i) =  aa(j,i) - bdi
        end do ! i
      end do ! i
c                       __                     _
c     Deviatoric term 2 mu [ I - 1/3 g x g ] / J

      mub1 = mub1 * trbb3
      mub2 = mub1 + mub1
      mub3 = mub2 * one3

      do i = 1,3
        aa(i  ,i  ) = aa(i  ,i  ) + mub2
        aa(i+3,i+3) = aa(i+3,i+3) + mub1
        do j = 1,3
          aa(i ,j ) = aa(i ,j )   - mub3
        end do ! i
      end do ! i

c     Tangent terms from second invariant

      nub2 = nub1 + nub1
      nub3 = two3 * nub2
      do i = 1,6
        bdi = nub3 * b2d(i)
        do j = 1,3
          aa(i,j) =  aa(i,j) - bdi
          aa(j,i) =  aa(j,i) - bdi
        end do ! j
      end do ! i

c     Last term: Form b in Cartesian tensor form

      b(1,1) = bb(1)
      b(1,2) = bb(4)
      b(1,3) = bb(6)
      b(2,1) = bb(4)
      b(2,2) = bb(2)
      b(2,3) = bb(5)
      b(3,1) = bb(6)
      b(3,2) = bb(5)
      b(3,3) = bb(3)

c     Term in tensor form
      dd = 0.0d0
      do jj = 1,9
        k = ia(jj)
        l = ja(jj)
        do ii = 1,9
          i = ia(ii)
          j = ja(ii)
          dd(ii,jj) = dd(ii,jj) + b(i,j)*b(k,l)
     &              - 0.5d0*(b(i,k)*b(j,l) + b(i,l)*b(j,k))
     &              + (one3*(de(i,k)*de(j,l) + de(i,l)*de(j,k))
     &              - four9* de(i,j)*de(k,l))*i2b
        end do ! ii
      end do ! jj

c     Compress to 6x6

      dd(:,4) = 0.5d0*(dd(:,4) + dd(:,5))
      dd(:,5) = 0.5d0*(dd(:,6) + dd(:,7))
      dd(:,6) = 0.5d0*(dd(:,8) + dd(:,9))

      dd(4,:) = 0.5d0*(dd(4,:) + dd(5,:))
      dd(5,:) = 0.5d0*(dd(6,:) + dd(7,:))
      dd(6,:) = 0.5d0*(dd(8,:) + dd(9,:))

      aa(1:6,1:6) = aa(1:6,1:6) + nub2*dd(1:6,1:6)

c     Compute inelastic parts for viscoelastic-damage model

      if(nint(d(40)).eq.2) then
        estore = d(22)*(1.d0 - d(23))*(i1b - 3.d0)*0.5d0
     &         + d(22)*d(23) *(i2b - 3.d0)*0.5d0   ! Used in damage
        call viscfd(d,f,finv,detf(1),hn(1),hn(2),hn(2+ntm),ntm, sig,aa)
      endif

c     Compute spatial deviatoric stresses and material moduli

      do i = 1,6
        sig(i) = sig(i) * detfi
        do j = 1,6
          aa(i,j) = aa(i,j) * detfi
        end do ! i
      end do ! i

c     Compute pressure and volumetric moduli

      jsw = nint(d(170))
      call fengy3(d,detf(3), u,up,upp,ha,hp,hpp,jsw)

c     Pressure and tangent (not mixed pressure)

      press =  up  + xlamd * hp
      upp   = (upp + xlamd * hpp) * detf(1)

c     Add volumetric correction to aa

      aa(1,1) = aa(1,1) - press + upp
      aa(1,2) = aa(1,2) + press + upp
      aa(1,3) = aa(1,3) + press + upp

      aa(2,1) = aa(2,1) + press + upp
      aa(2,2) = aa(2,2) - press + upp
      aa(2,3) = aa(2,3) + press + upp

      aa(3,1) = aa(3,1) + press + upp
      aa(3,2) = aa(3,2) + press + upp
      aa(3,3) = aa(3,3) - press + upp

      aa(4,4) = aa(4,4) - press
      aa(5,5) = aa(5,5) - press
      aa(6,6) = aa(6,6) - press

c     Check for structure vectors

      if(nint(d(260)).gt.0) then

        nstv = nint(d(260))
        do i = 1,nstv

          nprop = 268 + 3*nstv
          nvect = 258 + 3*nstv

c         Compute Holzapfel-Gasser model for fibers

          call pfiber(d(nprop),d(nvect),f(3,3), detf(1),sig, aa, 1)
        end do ! i

      endif

c     Add pressure terms

      do i = 1,3
        sig(i) = sig(i) + press
      end do ! i

c     Compute stored energy density

      engy  =  u + d(22)*(1.d0 - d(23))*(i1b - 3.d0)*0.5d0
     &                   + d(22)*d(23) *(i2b - 3.d0)*0.5d0

      end
