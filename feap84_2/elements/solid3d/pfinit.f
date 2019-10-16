c$Id:$
      subroutine pfinit(f,df,finv,detf, lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    27/04/2009
c       1. Initialize displacement gradient: dimension F=4. 10/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Initialize deformation gradient quantities

c      Inputs:
c        lint    - Number of quadrature points

c      Outputs:
c        f(9,4,*)  - Deformation gradient             = idenity
c        df(9,*)   - Incremental deformation gradient = zero
c        finv(9,*) - Inverse deformation gradient     = idenity
c        detf(4,*) - Determinant of deformation grad. = 1.0
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer    i,j,l,lint
      real*8     f(9,4,*), df(9,*), finv(9,*), detf(4,*)

      do l = 1,lint
        do i = 1,9
          do j = 1,4
            f(i,j,l)  = 0.0d0
          end do ! j
          finv(i,l) = 0.0d0
          df(i,l)   = 0.0d0
        end do ! i

        detf(1,l) = 1.d0
        detf(2,l) = 1.d0
        detf(3,l) = 0.d0
        detf(4,l) = 0.d0

        f(1,1,l)  = 1.0d0
        f(5,1,l)  = 1.0d0
        f(9,1,l)  = 1.0d0

        f(1,2,l)  = 1.0d0
        f(5,2,l)  = 1.0d0
        f(9,2,l)  = 1.0d0

        finv(1,l) = 1.0d0
        finv(5,l) = 1.0d0
        finv(9,l) = 1.0d0

      end do ! l

      end
