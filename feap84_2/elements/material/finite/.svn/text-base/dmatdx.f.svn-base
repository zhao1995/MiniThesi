c$Id:$
      subroutine dmatdx(dd,sigb,p_bar,p_mix)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute finite deformation mixed D arrays

c     Inputs:
c         dd(7,7)       Modified constitutive array
c         sigb(6)       Constitutive stresses
c         p_bar         Constitutive pressure
c         p_mix         Mixed pressure

c     Outputs:
c         dd(7,7)       Material tangent terms
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   i,j
      real*8    p_mix,p_bar,fac1,dd(7,7),sigb(6),sigb_d(6)

c     Compute deviator stress

      do i = 1,3
        sigb_d(i  ) = two3 *(sigb(i  ) - p_bar)
        sigb_d(i+3) = two3 * sigb(i+3)
      end do ! i

c     D_11: B_matrix part

      fac1 = p_mix - two3*p_bar
      do j = 1,3
        do i = 1,3
          dd(i,j) = dd(i,j) + fac1
        end do ! i
      end do ! j

      do j = 1,6
        do i = 1,3
          dd(i,j) = dd(i,j) - sigb_d(j)
          dd(j,i) = dd(j,i) - sigb_d(j)
        end do ! i
      end do ! j

      fac1 = p_bar - p_mix
      do j = 1,3
        dd(j  ,j  ) = dd(j  ,j  ) + fac1*2.d0
        dd(j+3,j+3) = dd(j+3,j+3) + fac1
      end do ! j

c     D_12: Coupling matrix with

      do j = 1,6
        dd(7,j) = dd(7,j) + sigb_d(j)
        dd(j,7) = dd(j,7) + sigb_d(j)
      end do ! j

c     D_22: Volumetric part

      dd(7,7) = dd(7,7) - one3*p_bar

      end
