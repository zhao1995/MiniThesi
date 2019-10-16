c$Id:$
      subroutine fbar1m(f,xji,theta,lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Form F-bar and left Cauchy-Green tensors

c     Inputs:
c        f(3,3)      - Deformation gradient
c        xji(2,*)    - Determinant of deformation gradient (J)
c        theta(2,*)  - Mixed determinant of deformation gradient
c        lint        - Number of quadrature points

c     Outputs:
c        f(3,3)      - Mixed deformation gradient
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   lint , i , l
      real*8    ration, ratio1, f(9,2,*), xji(2,*),theta(2,*)

c     Compute mixed deformation gradient

      do l = 1,lint
        ratio1 = (theta(1,l)/xji(1,l))**one3
        ration = (theta(2,l)/xji(2,l))**one3
        do i = 1,9
          f(i,1,l) = ratio1*f(i,1,l)
          f(i,2,l) = ration*f(i,2,l)
        end do ! i
      end do ! l

      end
