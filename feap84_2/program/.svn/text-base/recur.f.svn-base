c$Id:$
      subroutine recur (pn,dpn,pn1,x,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Recurrence relation for points

c      Inputs:
c         nn    - Order

c      Outputs:
c         pn    - Polynomial
c         dpn   - Derivative of polynomial
c         pn1   -
c         x     - Root estimat
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      integer    nn, n
      real*8     pn,dpn,pn1,x, c
      real*8     p1, p,dp, dp1, q,dq

      p1  = 1.d0
      p   = x
      dp1 = 0.d0
      dp  = 1.d0
      do n = 2,nn
        c   = 4.d0*dble(n-1)**4
     &       / (dble(2*n-1)*dble(2*n-3)*dble(2*n-2)**2)
        q   = x*p  - c*p1
        dq  = x*dp - c*dp1 + p
        p1  = p
        p   = q
        dp1 = dp
        dp  = dq
      end do ! n
      pn  = p
      dpn = dp
      pn1 = p1

      end
