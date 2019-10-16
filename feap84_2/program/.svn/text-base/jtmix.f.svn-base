c$Id:$
      subroutine jtmix(ea,eah,xa,xb,ra,rb,xbh,rah,rbh,
     &                 lambda,p,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form Rotational constraint arrays for connecting
c               rigid bodies.
c               Rotational constraint:  (xa - xb) dot ea = 0

c      Inputs:
c         ea(3)    - Orientation of axis at t_n+1
c         eah(3)   - Orientation of axis at t_n+alpha
c         xa(3)    - Location of body a connection point at t_n+1
c         xb(3)    - Location of body b connection point at t_n+1
c         ra(3)    - Location of center of mass of body 1 at t_n+1
c         rb(3)    - Location of center of mass of body 2 at t_n+1
c         xbh(3)   - Location of center of mass of body 1 at t_n+alpha
c         rah(3)   - Location of center of mass of body 2 at t_n+alpha
c         rbh(3)   - Location of center of mass of body 2 at t_n+alpha
c         lambda   - Lagrange multiplier value

c      Outputs:
c         p(*)     - Residual array for constraint
c         s(*,*)   - Tangent  array for constraint
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,nl
      real*8   p(3,6), s(17,17), lambda, lambdah
      real*8   xxedot, xrahedot, xrbehdot, xrbhedot
      real*8   ea(3), eah(3), xa(3), xb(3), xbh(3)
      real*8   ra(3), rah(3), rb(3), rbh(3)
      real*8   dx(3), xra(3), xrah(3), xrb(3), xrbh(3)
      real*8   dace(3), dbce(3), daceh(3), dbceh(3)
      real*8   tlama(3,3), tlamb(3,3), tlamc(3,3)

      save

      call pzero(p,18)
      call pzero(s,289)

      lambdah = 0.5d0*lambda

c     Form basic vectors:

c     dace  = (xb_n+1   - ra_n+1)   x ea_n+1
c     dbce  = (xb_n+1   - rb_n+1)   x ea_n+1
c     daceh = (xb_n+1/2 - ra_n+1/2) x ea_n+1/2
c     dbceh = (xb_n+1/2 - rb_n+1/2) x ea_n+1/2

c     tlama = lambda*(skew(xb_n+1/2 - ra_n+1/2)*skew(ea_n+1))/2
c     tlamb = lambda*(skew(ea_n+1/2)* skew(xb_n+1  - rb_n+1))/2
c     tlamc = lambda*(skew(xb_n+1/2 - rb_n+1/2)*skew(ea_n+1))/2

c     dx   = xa_n+1   - xb_n+1
c     xra  = xb_n+1   - ra_n+1
c     xrah = xb_n+1/2 - ra_n+1/2
c     xrb  = xb_n+1   - rb_n+1
c     xrbh = xb_n+1/2 - rb_n+1/2

c     xxedot   = (xa_n+1 - xb_n+1)     * ea_n+1
c     xrahedot = (xb_n+1/2 - ra_n+1/2) * ea_n+1
c     xrbehdot = (xb_n+1 - rb_n+1)     * ea_n+1/2
c     xrbhedot = (xb_n+1/2 - rb_n+1/2) * ea_n+1

      do i = 1,3
        dx(i)   = xa(i)  - xb(i)
        xra(i)  = xb(i)  - ra(i)
        xrah(i) = xbh(i) - rah(i)
        xrb(i)  = xb(i)  - rb(i)
        xrbh(i) = xbh(i) - rbh(i)
      end do ! i

      daceh(1) = xrah(2)*eah(3) - xrah(3)*eah(2)
      daceh(2) = xrah(3)*eah(1) - xrah(1)*eah(3)
      daceh(3) = xrah(1)*eah(2) - xrah(2)*eah(1)

      dbceh(1) = xrbh(2)*eah(3) - xrbh(3)*eah(2)
      dbceh(2) = xrbh(3)*eah(1) - xrbh(1)*eah(3)
      dbceh(3) = xrbh(1)*eah(2) - xrbh(2)*eah(1)

      dace(1) = xra(2)*ea(3) - xra(3)*ea(2)
      dace(2) = xra(3)*ea(1) - xra(1)*ea(3)
      dace(3) = xra(1)*ea(2) - xra(2)*ea(1)

      dbce(1) = xrb(2)*ea(3) - xrb(3)*ea(2)
      dbce(2) = xrb(3)*ea(1) - xrb(1)*ea(3)
      dbce(3) = xrb(1)*ea(2) - xrb(2)*ea(1)

      xxedot   = dx(1)*ea(1)   + dx(2)*ea(2)   + dx(3)*ea(3)
      xrahedot = xrah(1)*ea(1) + xrah(2)*ea(2) + xrah(3)*ea(3)
      xrbhedot = xrbh(1)*ea(1) + xrbh(2)*ea(2) + xrbh(3)*ea(3)
      xrbehdot = xrb(1)*eah(1) + xrb(2)*eah(2) + xrb(3)*eah(3)

      do i = 1,3
        do j = 1,3
          tlama(i,j) = xrah(j)*ea(i)
          tlamb(i,j) = eah(j) *xrb(i)
          tlamc(i,j) = xrbh(j)*ea(i)
        end do ! j
        tlama(i,i) = tlama(i,i) - xrahedot
        tlamb(i,i) = tlamb(i,i) - xrbehdot
        tlamc(i,i) = tlamc(i,i) - xrbhedot
      end do ! i

c     Form residual

      do i = 1,3
        p(i,1) = -lambda*eah(i)
        p(i,2) = -lambda*daceh(i)
        p(i,3) =  lambda*eah(i)
        p(i,4) =  lambda*dbceh(i)
      end do ! i
      p(1,5) = -xxedot

c     Form Tangent

      s(1,5) =  lambdah*ea(3)
      s(1,6) = -lambdah*ea(2)
      s(2,4) = -lambdah*ea(3)
      s(3,4) =  lambdah*ea(2)
      s(2,6) =  lambdah*ea(1)
      s(3,5) = -lambdah*ea(1)

      s(4,2) = -lambdah*eah(3)
      s(4,3) =  lambdah*eah(2)
      s(5,1) =  lambdah*eah(3)
      s(6,1) = -lambdah*eah(2)
      s(5,3) = -lambdah*eah(1)
      s(6,2) =  lambdah*eah(1)

      nl = 13
      do i = 1,3

        s(i  ,nl) =  eah(i)
        s(i+3,nl) =  daceh(i)
        s(i+6,nl) = -eah(i)
        s(i+9,nl) = -dbceh(i)
        s(nl,i  ) =  ea(i)
        s(nl,i+3) =  dace(i)
        s(nl,i+6) = -ea(i)
        s(nl,i+9) = -dbce(i)

        do j = 1,3
          s(i+6,j+3) = -s(i,j+3)
          s(i+3,j+6) = -s(i+3,j)

          s(i+3,j+3) = -tlama(i,j)*lambdah
          s(i+3,j+9) =  tlamb(i,j)*lambdah
          s(i+9,j+3) =  tlamc(i,j)*lambdah
          s(i+9,j+9) = -tlamb(i,j)*lambdah
        end do ! j

      end do ! i

      end
