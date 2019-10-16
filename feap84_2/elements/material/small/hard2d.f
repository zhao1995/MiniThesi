c$Id:$
      function hard2d(d,ep,xkprim,hprim)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Linear-saturation hardeing model for plasticity

c      Inputs:
c         d(*)     - Material parameters
c         ep       - Accumlated plastic strain value

c      Outputs:
c         xkprim   - Kinematic hardening part
c         hprim    - Isotropic hardentin part
c         hard2d   - Value of hardening parameter
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      real*8    hard2d,y0,yinf,ylin,beta,ep,eep,xk,xkprim,hprim,d(*)

c     Extract constants

      y0     = d(41)
      yinf   = d(42)
      beta   = d(43)
      ylin   = d(44)
      hprim  = d(45)

c     Compute yield function and derivative

      eep    = exp(-beta*ep)
      xk     = y0 + (yinf - y0)*(1.d0 - eep) + ylin*ep
      xkprim = beta*(yinf - y0)*eep          + ylin

      hard2d = xk*xk*one3

      end
