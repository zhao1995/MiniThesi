c$Id:$
      subroutine arruda(d,detf,bb, sig,dd, xlamd,ha, estore)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'pconstant.h'                 14/11/2006
c       2. Dimension detf(3), use in fengy                  22/11/2011
c          Modify to use bb - 1.0
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compressible Arruda-Boyce model
c                (J_2/3 regularization)
c                     __                   _
c                 W(J,be) = U(J) + 0.5*mu*[I_1 - 3)
c                                        _   _
c                                + m/10*(I_1*I_1 - 9)
c                                              _   _   _
c                                + 11*m*m/525*(I_1*I_1*I_1 - 27)]
c                 _                  __
c                 I_1     = J^(-2/3)(be:1 - 3)

c     Input:
c          d(*)    -  Material parameters
c          detf(3) -  Determinant of deforamtion gradient
c          bb(6)   -  Left Cauchy deformation tensor
c          xlamd   -  Augmented "penalty" value

c     Output:
c          sig(*)  -  Stresses at point.
c                     N.B. 1-d models use only sig(1)
c          dd(6,*) -  Current material tangent moduli
c                     N.B. 1-d models use only dd(1,1) and dd(2,1)
c          ha      -  Augmented function
c          estore  -  Stored energy value
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include  'pconstant.h'
      include  'elcount.h'

      integer   i,j, jsw
      real*8    bb(6),be(6), detf(3), d(*), sig(*),dd(6,*)

      real*8    detfi, j23, trbe3, bei, mub1,mub2,mub3,mub4, estore
      real*8    u, up, upp, ha, hp, hpp, mm, press, xlamd

      save

c     Compute deviatoric be

      detfi = 1.d0/detf(1)
      j23   = detfi**two3
      do i = 1,6
        be(i) = bb(i) * j23
      end do ! i

      trbe3 = (be(1) + be(2) + be(3)) * one3
      be(1) = be(1) - trbe3
      be(2) = be(2) - trbe3
      be(3) = be(3) - trbe3
      trbe3 = trbe3 + j23

c     Compute deviatoric Kirchhoff stress tensor.

      mub1 = d(22) ! Arruda-Boyce mu
      mm   = d(23) ! Arruda-Boyce m = 1/n
      mub4 = mub1*mm*(0.4d0 + 132.d0*mm*trbe3/175.d0)
      mub1 = mub1*(1.d0 + 0.6d0*mm*trbe3*(1.d0 + 33.d0/35.d0*mm*trbe3))
      do i = 1,6
        sig(i) = mub1 * be(i)
      end do ! i

c     Compute tangent tensor

      mub3 = two3 * mub1
      do i = 1,6
        bei = mub4 * be(i)
        do j = 1,6
          dd(i,j) = bei*be(j)
        end do ! j
      end do ! i
      do i = 1,6
        bei = mub3 * be(i)
        do j = 1,3
          dd(i,j) =  dd(i,j) - bei
          dd(j,i) =  dd(j,i) - bei
        end do ! j
      end do ! i

      mub1 = mub1 * trbe3
      mub2 = mub1 + mub1
      mub3 = mub2 * one3

      do i = 1,3
        dd(i  ,i  ) = dd(i  ,i  ) + mub2
        dd(i+3,i+3) = dd(i+3,i+3) + mub1
        do j = 1,3
          dd(i ,j ) = dd(i ,j )   - mub3
        end do ! j
      end do ! i

c     Compute spatial deviatoric stresses and material moduli

      do i = 1,6
        sig(i) = sig(i) * detfi
        do j = 1,6
          dd(i,j) = dd(i,j) * detfi
        end do ! j
      end do ! i

c     Compute pressure and volumetric moduli

      jsw = nint(d(170))
      call fengy3(d,detf(3), u,up,upp,ha,hp,hpp,jsw)

c     Pressure and its tangent

      press =   up  + xlamd * hp
      upp   =  (upp + xlamd * hpp) * detf(1)

c     Add volumetric correction to dd

      dd(1,1) = dd(1,1) - press + upp
      dd(1,2) = dd(1,2) + press + upp
      dd(1,3) = dd(1,3) + press + upp

      dd(2,1) = dd(2,1) + press + upp
      dd(2,2) = dd(2,2) - press + upp
      dd(2,3) = dd(2,3) + press + upp

      dd(3,1) = dd(3,1) + press + upp
      dd(3,2) = dd(3,2) + press + upp
      dd(3,3) = dd(3,3) - press + upp

      dd(4,4) = dd(4,4) - press
      dd(5,5) = dd(5,5) - press
      dd(6,6) = dd(6,6) - press

c     Add pressure terms

      do i = 1,3
        sig(i) = sig(i) + press
      end do ! i

c     Energy

      estore = 0.0d0  !  NEED TO ADD

      end
