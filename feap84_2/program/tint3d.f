c$Id:$
      subroutine tint3d(ll,lint,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constant values from 'pconstant.h'           07/02/2007
c          Correct weight on 5-pt formula (ll = 3)
c       2. Make sum of weights 1.0                          13/08/2007
c       3. Add 14 point formula                             24/09/2007
c       4. Add  8 point formula; compute 14 point values    29/12/2007
c       5. Add 10 and 15 point formulas                     08/01/2008
c       6. Correct point for 11 point (ll=4) formula        21/03/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:   Gauss quadrature for 3-d tetrahedral element

c      Reference: M. Gellert and R. Harbord, "Moderate degree cubature
c                 formulas for 3-D tetrahedral finite-element approx-
c                 imation", Com. in Appl. Num. Meth., Vol. 7, 1991,
c                 pp 487-495.

c      Inputs:
c         ll       - Type of quadrature

c      Outputs:
c         lint     - Number of quadrature points
c         s(5,*)   - Values of volume coordinates and weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   i, j, ll, lint
      real*8    s(5,*), alp(2), con

      integer  ind1(6),ind2(6),ind3(6),ind4(6)

      save

      data     ind1/ 1,2,2,1,1,2/
      data     ind2/ 1,1,2,2,2,1/
      data     ind3/ 2,1,1,2,1,2/
      data     ind4/ 2,2,1,1,2,1/

c     1 pt. quadrature O(h^2)

      if(ll.eq.1) then
        lint = 1
        do i = 1,4
          s(i,1) = 0.25d0
        end do ! i
        s(5,1) = 1.0d0

c     4 pt. quadrature O(h^2) - nodes of linear element

      elseif(ll.eq.-1) then
        lint = 4
        s(5,4) = 0.25d0
        do i = 1,4
          do j = 1,4
            s(i,j) = 0.0d0
          end do ! j
          s(i,i) = 1.d0
          s(5,i) = s(5,4)
        end do ! i

c     4 pt. quadrature O(h^3)

      elseif(ll.eq.2) then
        lint = 4

        con    = sqrt(5.d0)
        alp(1) = 0.25d0 + 0.15d0*con
        alp(2) = 0.25d0 - 0.05d0*con

        s(5,4) = 0.25d0
        do i = 1,4
          do j = 1,4
            s(i,j) = alp(2)
          end do ! j
          s(i,i) = alp(1)
          s(5,i) = s(5,4)
        end do ! i

c     5  pt. quadrature O(h^4) -- has negative weight

      elseif(ll.eq.3) then
        lint = 5
        do i = 1,4
          s(i,1) = one6
          s(i,2) = s(i,1)
          s(i,3) = s(i,1)
          s(i,4) = s(i,1)
          s(i,5) = 0.25d0
        end do ! i
        do i = 1,4
          s(i,i) = 0.5d0
        end do ! i
        s(5,1) =  0.45d0
        s(5,2) =  s(5,1)
        s(5,3) =  s(5,1)
        s(5,4) =  s(5,1)
        s(5,5) = -0.80d0

c     11 pt. quadrature O(h^4) -- has negative weight

      elseif(ll.eq.4) then
        lint = 11

        con    = 1.d0/14.d0
        alp(1) = 0.25d0 + 0.25d0*sqrt(5.d0*con)
        alp(2) = 0.50d0 - alp(1)

        do i = 1,4
          do j = 1,4
            s(j,i) = con
          end do ! j
          s(i,i) = 11.d0*con
        end do ! i

        con = 343.d0/7500.d0
        do i = 1,4
          s(5,i) = con
        end do ! i

        con = 56.d0/375.d0
        do i = 1,6
          s(1,i+4) = alp(ind1(i))
          s(2,i+4) = alp(ind2(i))
          s(3,i+4) = alp(ind3(i))
          s(4,i+4) = alp(ind4(i))
          s(5,i+4) = con
        end do ! i
        s(1,11) =  0.25d+00
        s(2,11) =  0.25d+00
        s(3,11) =  0.25d+00
        s(4,11) =  0.25d+00
        s(5,11) = -148.d0/1875.d0

c     11 pt. quadrature O(h^4) -- has no negative weight

      elseif(ll.eq.-4) then
        lint = 11
        do i = 1,4
          do j = 1,10
            s(i,j) = 0.0d0
          end do ! j
          s(i, i  ) = 1.00d0
          s(i, i+4) = 0.50d0
          s(i, i+7) = 0.50d0
          s(i, 11 ) = 0.25d0
        end do ! i
        s(2, 5) = 0.50d0
        s(3, 6) = 0.50d0
        s(1, 7) = 0.50d0
        do j = 1,4
          s(5,j) = 1.d0/60.d0
        end do ! j
        do j = 5,10
          s(5,j) = 1.d0/15.d0
        end do ! j
        s(5,11) = 24.d0/45.d0

c     8 pt. quadrature O(h^3)

      elseif(ll.eq. 8) then

        lint = 8
        do i = 1,4
          do j = 1,3
            s(j,i  ) = 0.0d0
            s(j,i+4) = one3
          end do ! j
          s(i,i  ) = 1.0d0
          s(i,i+4) = 0.0d0
          s(5,i  ) = 0.025d0
          s(5,i+4) = 0.225d0
        end do ! i

c     10 pt. quadrature O(h^4)

      elseif(ll.eq.10) then

        lint = 10

        do i = 1,4
          do j = 1,4
            s(j,i) = b10
          end do ! j
          s(i,i) = a10
          s(5,i) = p10
        end do ! i

        s(1, 5) = 0.5d0
        s(2, 5) = 0.5d0
        s(3, 5) = 0.0d0

        s(1, 6) = 0.0d0
        s(2, 6) = 0.5d0
        s(3, 6) = 0.5d0

        s(1, 7) = 0.5d0
        s(2, 7) = 0.0d0
        s(3, 7) = 0.5d0

        s(1, 8) = 0.5d0
        s(2, 8) = 0.0d0
        s(3, 8) = 0.0d0

        s(1, 9) = 0.0d0
        s(2, 9) = 0.5d0
        s(3, 9) = 0.0d0

        s(1,10) = 0.0d0
        s(2,10) = 0.0d0
        s(3,10) = 0.5d0

        do i = 5,10
          s(5,i) = q10
        end do ! i

c     14 pt. quadrature O(h^5)

      elseif(ll.eq.14 .or. ll.eq. 5) then

        lint = 14

        do i = 1,4
          do j = 1,3
            s(j,i  ) = b14
            s(j,i+4) = d14
          end do ! j
          s(i,i  ) = a14
          s(i,i+4) = c14
          s(5,i  ) = p14
          s(5,i+4) = q14
        end do ! i

        do i = 9,14
          do j = 1,3
            s(j,i) = f14
          end do ! j
          s(5,i) = r14
        end do ! i
        s(1, 9) = e14
        s(2, 9) = e14
        s(2,10) = e14
        s(3,10) = e14
        s(3,11) = e14
        s(1,12) = e14
        s(3,12) = e14
        s(2,13) = e14
        s(1,14) = e14

c     15 pt. quadrature O(h^6)

      elseif(ll.eq.15) then

        lint = 15

        do i = 1,4
          do j = 1,3
            s(j,i  ) = 1.0d0/11.0d0
            s(j,i+4) = one3
          end do ! j
          s(i,i  ) = 8.0d0/11.0d0
          s(i,i+4) = 0.0d0
          s(5,i  ) = 0.116452490860289742d-01
          s(5,i+4) = 0.602678571428571597d-02
        end do ! i

        do i = 9,14
          do j = 1,3
            s(j,i) = 0.433449846426335728d+00
          end do ! j
          s(5,i) = 0.109491415613864534d-01
        end do ! i
        s(3, 9) = 0.66550153573664281d-01
        s(1,10) = 0.66550153573664281d-01
        s(2,11) = 0.66550153573664281d-01
        s(2,12) = 0.66550153573664281d-01
        s(3,12) = 0.66550153573664281d-01
        s(1,13) = 0.66550153573664281d-01
        s(3,13) = 0.66550153573664281d-01
        s(1,14) = 0.66550153573664281d-01
        s(2,14) = 0.66550153573664281d-01

        do j = 1,3
          s(j,15) = 0.25d0
        end do ! j
        s(5,15) = 0.302836780970891856d-01

c     16 pt. quadrature O(h^5)

      else
        lint = 16
        s(5,4) = 0.8395632516687135d-02*6.d0
        do i = 1,3
          do j = 1,4
            s(i,j) = 0.7611903264425430d-01
          end do ! j
          s(i,i) = 0.7716429020672371d+00
          s(5,i) = s(5,4)
        end do ! i
        do i = 5,16
          s(5,i) = 0.1109034477221540d-01*6.d0
        end do ! i

        s(1, 5) = 0.1197005277978019d+00
        s(2, 5) = 0.7183164526766925d-01
        s(3, 5) = 0.4042339134672644d+00
        s(1, 6) = 0.4042339134672644d+00
        s(2, 6) = 0.1197005277978019d+00
        s(3, 6) = 0.7183164526766925d-01
        s(1, 7) = 0.4042339134672644d+00
        s(2, 7) = 0.4042339134672644d+00
        s(3, 7) = 0.1197005277978019d+00
        s(1, 8) = 0.7183164526766925d-01
        s(2, 8) = 0.4042339134672644d+00
        s(3, 8) = 0.4042339134672644d+00

        s(1, 9) = 0.1197005277978019d+00
        s(2, 9) = 0.4042339134672644d+00
        s(3, 9) = 0.7183164526766925d-01

        s(1,10) = 0.4042339134672644d+00
        s(2,10) = 0.1197005277978019d+00
        s(3,10) = 0.4042339134672644d+00

        s(1,11) = 0.7183164526766925d-01
        s(2,11) = 0.4042339134672644d+00
        s(3,11) = 0.1197005277978019d+00

        s(1,12) = 0.4042339134672644d+00
        s(2,12) = 0.7183164526766925d-01
        s(3,12) = 0.4042339134672644d+00

        s(1,13) = 0.1197005277978019d+00
        s(2,13) = 0.4042339134672644d+00
        s(3,13) = 0.4042339134672644d+00

        s(1,14) = 0.7183164526766925d-01
        s(2,14) = 0.1197005277978019d+00
        s(3,14) = 0.4042339134672644d+00

        s(1,15) = 0.4042339134672644d+00
        s(2,15) = 0.7183164526766925d-01
        s(3,15) = 0.1197005277978019d+00

        s(1,16) = 0.4042339134672644d+00
        s(2,16) = 0.4042339134672644d+00
        s(3,16) = 0.7183164526766925d-01

      endif

c     Compute fourth points

      do j = 1,lint
        s(4,j) = 1.d0 - (s(1,j) + s(2,j) + s(3,j))
      end do ! j

      end
