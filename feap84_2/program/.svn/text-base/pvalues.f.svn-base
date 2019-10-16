c$Id:$
      subroutine pvalues()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    14/11/2006
c       1. Add 'one12'                                      09/03/2007
c       2. Add 'one36'                                      29/03/2007
c       3. Add 14-pt quadrature values                      29/12/2007
c       4. Add 10-pt quadrature values                      08/01/2008
c       5. Add pi180 for pi/180                             22/02/2009
c       6. Add one4,one5, fac1r to fac6r                    08/12/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute values of constants used in program

c      Inputs:
c         none

c      Outputs:
c         Values are output through common 'pcommon.h'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      real*8    c1,c2,c3,c4

      one2  = 0.5d0                   ! One  half

      one3  = 1.0d0
      one3  = one3/3.0d0              ! One  third
      two3  = one3*2.0d0              ! Two  thirds
      four3 = two3*2.0d0              ! Four thirds

      one4  = 0.25d0                  ! One  fourth
      one5  = 0.20d0                  ! One  fifthe

      one6  = one3*0.5d0              ! One  sixth
      one12 = one6*0.5d0              ! One  twelfth
      one36 = one12*one3              ! One  thirty-sixth

      one9  = one3**2                 ! One  ninth
      four9 = one9*4.0d0              ! Four ninths

      pi    = acos(-1.d0)             ! pi
      pi23  = (pi*2.0d0)/3.0d0        ! pi*2/3
      pi180 = pi/180.0d0              ! pi/180

      sqrt2 = 2.0d0
      sqrt2 = sqrt(sqrt2)             ! Square root of 2

      sqt13 = 1.0d0
      sqt13 = sqt13/3.0d0
      sqt13 = sqrt(sqt13)             ! Square root of 1/3

      sqt23 = 2.0d0
      sqt23 = sqt23/3.0d0
      sqt23 = sqrt(sqt23)             ! Square root of 2/3

      sqtp6 = 0.6d0
      sqtp6 = sqrt(sqtp6)             ! Sqare root of 0.6

      sqt48 = 4.8d0
      sqt48 = sqrt(sqt48)             ! Sqare root of 4.8

      five9 = 5.0d0
      five9 = five9/9.0d0             ! 5/9

      eight9 = 8.0d0
      eight9 = eight9/9.0d0           ! 8/9

      thty29 = 32.0d0
      thty29 = thty29/9.0d0           ! 32/9

c     Parameters for 14-pt quadrature on tetrahedra

      c1 = 1.d0/(46.d0*sqrt(46.d0))
      c2 = acos(c1) + two3*asin(c1)
      c3 = (104.d0 + 8.d0*sqrt(46.d0)*cos(c2))*one3
      c4 = sqrt(49.d0 - c3)

      b14 = (7.d0 + c4)/c3
      a14 = 1.d0 - 3.d0*b14
      d14 = (7.d0 - c4)/c3
      c14 = (1.d0 - 3.d0*d14)
      p14 = (98.d0 - c3 - 14.d0*c4)/(1680.d0*c4*(b14 - a14)**3)
      q14 = (98.d0 - c3 + 14.d0*c4)/(1680.d0*c4*(c14 - d14)**3)
      r14 = (1.d0 - 4.d0*(p14+q14))*one6
      e14 = (1.d0 + (2.d0/(105.d0*r14))**0.25d0)*0.25d0
      f14 = (1.d0 - 2.d0*e14)*one2

c     Parameters for 10-pt quadrature on tetrahedra

      a10 = (1.d0 + (17.d0 + 12.d0*sqrt2)**one3
     &            + (17.d0 - 12.d0*sqrt2)**one3)*0.125d0
      b10 = (1.d0 - a10)*one3
      p10 = 1.d0/(120.d0*b10*(3.d0 - 8.d0*b10))
      p10 = 0.2177650698804054d0
      q10 = (1.d0 - 4.d0*p10)*one6

c     Reciprocal factorials

      fac1r = 1.0d0
      fac2r = 0.5d0
      fac3r = one6
      fac4r = fac3r*one4
      fac5r = fac4r*one5
      fac6r = fac5r*one6

      end
