c$Id:$
      subroutine pstres(sig,p1,p2,p3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute principal stresses for 2-d problems.
c               N.B. Should replace with 'pstr2d' when possible

c      Input:
c         sig(1) - Stresses in order: sig-xx, sig-xy, sig-yy
c      Output:
c         p1     - Principal stresse sig-1
c         p2     - Principal stresse sig-2
c         p3     - angle (degrees): sig-xx to sig-1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      real*8    sig(3), xi1, xi2, rho, p1, p2, p3

      save

      xi1  = (sig(1) + sig(3))*0.5d0
      xi2  = (sig(1) - sig(3))*0.5d0
      rho  = sqrt(xi2*xi2 + sig(2)*sig(2))
      p1   = xi1 + rho
      p2   = xi1 - rho
      if(xi2.ne.0.0d0) then
        p3 = 90.0d0*atan2(sig(2),xi2)/pi
      else
        p3 = 45.0d0
      endif

      end
