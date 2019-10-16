c$Id:$
      real*8 function yldgpl(ep,salm,saltrm,g,yield,beta,delta,
     &                       hiso,hkin,tt,lambda)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Generalized plasticity yield

c     Inputs:
c          ep     - Value of plastic strain at time t(n)
c          salm   - Norm of deviator less back stress at time t(n)
c          saltrm - Norm of trial deviator less back stress
c          g      - Shear modulus
c          yield  - Initial yield
c          beta   - Increase to limit yield * sqrt(2/3)
c          delta  - Growth rate * 2/3
c          hiso   - Isotropic hardening
c          hkin   - Kinematic hardening
c          tt     - Sqrt(2/3)

c     Outputs:
c          yldgpl - yield function value for current iterate
c          lambda - consistency parameter (from constitutive model)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pconstant.h'

      real*8    salm,saltrm,ep,g,yield,beta,delta,hiso,hkin,tt, lambda
      real*8    a1,a2,a3,aa,bb,discr,g2

      save

c     Set parameters

      g2  = g + (hiso + hkin)*one3

c     Solution of quadratic equation

      a1 = saltrm - tt * (yield + hiso*ep)
      a2 = saltrm - salm
      a3 = 2.0d0 * g - delta

      aa =   2.0d0 * g2 * a3
      bb = - a1 * a3 - 2.0d0 * a2 * g2
     &     - (delta + two3*(hkin + hiso))*beta

      discr = bb * bb - 4.0d0 * aa * a1 * a2

      if (discr.lt.0.0d0) then
        write(iow,*) ' DISCR ERROR ', discr
      end if

      lambda = -0.5d0 * (bb + sqrt(abs(discr)))/aa
      yldgpl =     tt * (yield + hiso*(ep + tt*lambda))

      end
