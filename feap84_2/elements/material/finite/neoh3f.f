c$Id:$
      subroutine neoh3f(d,f,finv,detf,bb, hn,ntm,sig,aa,
     &                  xlamd,ha,engy)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'pconstant.h'                 14/11/2006
c       2. Cast d(40) as integer for compare                07/06/2007
c       3. Convert to use displacement gradient             21/11/2011
c       4. Dimension detf(3), use in fengy                  22/11/2011
c       5. Pass detf(1) to viscfd                           01/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compressible Neohookean model
c                (J_2/3 regularization)
c                   _ __      _            __
c                 W(J,be) = U(J) + 0.5*mu*(be:1 - 3)
c                     __
c                     be  = J^{2/3} * be

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'elengy.h'
      include  'pconstant.h'

      integer   i, j, ntm, jsw
      real*8    detf(3), detfi, j23, trbe3, bei, mub1, mub2, mub3
      real*8    u, up, upp, ha, hp, hpp, press, xlamd, engy
      real*8    d(*), f(3,3), finv(3,3), bb(6), hn(*), sig(6), aa(6,6)
      real*8    be(6)

      save

c     Compute deviatoric be

      detfi = 1.d0/detf(1)
      j23   = detfi**two3
      do i = 1,ntm
        be(i) = bb(i) * j23
      end do ! i

      trbe3  = (be(1) + be(2) + be(3)) * one3
      be(1) = be(1) - trbe3
      be(2) = be(2) - trbe3
      be(3) = be(3) - trbe3

c     Compute deviatoric Kirchhoff stress tensor.

      mub1 = d(22)
      do i = 1,ntm
        sig(i) = mub1 * be(i)
      end do ! i

c     Compute tangent tensor
c                                  __             __     _
c     Rank one update: -2/3 mu * ( be x g +  g x  be ) / J

      mub3 = two3 * mub1
      do i = 1,ntm
        bei = mub3 * be(i)
        do j = 1,3
          aa(i,j) =  aa(i,j) - bei
          aa(j,i) =  aa(j,i) - bei
        end do ! j
      end do ! i
c                       __                     _
c     Deviatoric term 2 mu [ I - 1/3 g x g ] / J

      mub1 = mub1 *(trbe3 + j23)
      mub2 = mub1 + mub1
      mub3 = mub2 * one3

      do i = 1,3
        aa(i  ,i  ) = aa(i  ,i  ) + mub2
        aa(i+3,i+3) = aa(i+3,i+3) + mub1
        do j = 1,3
          aa(i ,j ) = aa(i ,j )   - mub3
        end do ! j
      end do ! i

c     Compute inelastic parts

c     1. Plasticity

      if(nint(d(40)).eq.1) then

c     2. Viscoelasticity

      elseif(nint(d(40)).eq.2) then

        estore = d(22)*trbe3*1.5d0  ! Used for damage model
        call viscfd(d,f,finv,detf(1),hn(1),hn(2),hn(2+ntm),ntm, sig,aa)

      endif

c     Compute spatial deviatoric stresses and material moduli

      do i = 1,ntm
        sig(i) = sig(i) * detfi
        do j = 1,ntm
          aa(i,j) = aa(i,j) * detfi
        end do ! j
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

c     Add pressure terms

      do i = 1,3
        sig(i) = sig(i) + press
      end do ! i

c     Compute stored energy density

      engy = u + d(22)*trbe3*1.5d0

      end
