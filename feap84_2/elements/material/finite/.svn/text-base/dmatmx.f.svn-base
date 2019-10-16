c$Id:$
      subroutine dmatmx ( aa, dd )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Project 6x6 AA-matrix onto 7x7 DD-matrix for mixed model
c              (N.B. Jacobian ratio is in dvol)
c                     | D_11   D_12 |
c                DD = |             |
c                     | D_21   D_22 |
c              where:
c                D_11  = 6 x 6 Deviatoric part of matrix
c                D_12  = 6 x 1 Coupling   part of matrix
c                D_21  = 1 x 6 Coupling   part of matrix
c                D_22  = 1 x 1 Volumetric part of matrix

c     Inputs:
c        aa(6,6)   - Material tangent matrix (based on F)

c     Outputs:
c        dd(7,7)   - Mixed material tangent for stiffness computations.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   i , j
      real*8    aa(6,6), dd(7,7)

c     Load moduli from constitution

      do i = 1,6
        do j = 1,6
          dd(i,j) = aa(i,j)
        end do ! j
      end do ! i

c     Compute left and right multiples with trace

      do i = 1,6
        dd(i,7) = (aa(i,1) + aa(i,2) + aa(i,3))*one3
        dd(7,i) = (aa(1,i) + aa(2,i) + aa(3,i))*one3
      end do ! i

c     Convert upper 6 x 6 to a deviatoric D_11

      do i = 1,6
        do j = 1,3
          dd(i,j) = dd(i,j) - dd(i,7)
          dd(j,i) = dd(j,i) - dd(7,i)
        end do ! j
      end do ! i

c     Form last term, D_22

      dd(7,7) = (dd(1,7) + dd(2,7) + dd(3,7))*one3

c     Final update to form D_11, D_12 and D_21

      do i = 1,3
        dd(i,7) = dd(i,7) - dd(7,7)
        dd(7,i) = dd(7,i) - dd(7,7)
        do j = 1,3
          dd(i,j) = dd(i,j) + dd(7,7)
        end do ! j
      end do ! i

      end
