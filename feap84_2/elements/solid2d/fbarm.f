c$Id:$
      subroutine fbarm(f,detf,theta,lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Get 'one3' from 'pconstant.h'                    16/09/2009
c       2. Dimension f, detf and theta to 4.                21/11/2011
c          Add computation of G-bar
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Form F-bar and G-bar

c     Inputs:
c        f(3,3)      - Deformation gradient
c        detf(4,*)   - Determinant of deformation gradient (J)
c        theta(4,*)  - Mixed determinant of deformation gradient
c        lint        - Number of quadrature points

c     Outputs:
c        f(9,4,*)      - Mixed deformation gradient: Voigt form
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   lint , i , l
      real*8    ration, ratio1, psin_d, psi1_d, psin_m, psi1_m
      real*8    f(9,4,*), detf(4,*),theta(4,*)

c     External  function

      real*8    ddetj

c     Loop over quadrature points

      do l = 1,lint

c       Compute mixed deformation gradient: F-bar

        psi1_d = abs(  detf(1,l))**one3
        ratio1 = abs(theta(1,l))**one3/psi1_d

        psin_d = abs(  detf(2,l))**one3
        ration = abs(theta(2,l))**one3/psin_d

        do i = 1,9
          f(i,1,l) = ratio1*f(i,1,l)
          f(i,2,l) = ration*f(i,2,l)
        end do ! i

c       Compute mixed displacement gradient: G-bar

        if(abs(theta(3,l) - detf(3,l)) .lt. 0.001d0) then
          psi1_m = ddetj(detf(3,l),theta(3,l))/psi1_d
        else
          psi1_m = ratio1 - 1.0d0
        endif
        if(abs(theta(4,l) - detf(4,l)) .lt. 0.001d0) then
          psin_m = ddetj(detf(4,l),theta(4,l))/psin_d
        else
          psin_m = ration - 1.0d0
        endif
        do i = 1,9
          f(i,3,l) = f(i,3,l) + psi1_m*f(i,3,l)
          f(i,4,l) = f(i,4,l) + psin_m*f(i,4,l)
        end do ! i

c       Add effect of identity to G-bar

        f(1,3,l) = f(1,3,l) + psi1_m
        f(5,3,l) = f(5,3,l) + psi1_m
        f(9,3,l) = f(9,3,l) + psi1_m

        f(1,4,l) = f(1,4,l) + psin_m
        f(5,4,l) = f(5,4,l) + psin_m
        f(9,4,l) = f(9,4,l) + psin_m
      end do ! l

      end
