c$Id:$
      subroutine legendre1d(xi, p, leg, dleg)
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute value of order p Legendre polynomial at xi

c      Recursion:
c        (n+1)* P_{n+1} = (2*n+1)*xi*P_n          - n* P_{n-1}
c        (n+1)*dP_{n+1} = (2*n+1)*(P_n + xi*dP_n) - n*dP_{n-1}

c      Inputs:
c        xi     - Point to evaluate
c        p      - Order of legendre polynomial

c      Outputs:
c        leg    - Legendre polynomial
c        dleg   - Derivative of legendre polynomial
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    p, n
      real*8     xi, leg, dleg, ln, lm,  dn, dm

c     Case: p = 0

      if(p.eq.0) then
        leg =  1.0d0
        dleg = 0.0d0

c     Case: p = 1

      elseif(p.eq.1) then
        leg =  xi
        dleg = 1.0d0

c     Perform recursion for order p

      else
        lm = 1.0d0
        ln = xi
        dm = 0.0d0
        dn = 1.0d0
        do n = 2,p
          leg  = (dble(2*n+1)*xi*ln - dble(n)*lm)/dble(n+1)
          dleg = (dble(2*n+1)*(ln + xi*dn) - dble(n)*dm)/dble(n+1)
          lm   = ln
          dm   = dn
          ln   = leg
          dn   = dleg
        end do ! n
      endif

      end
