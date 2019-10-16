c$Id:$
      subroutine quatex(t, q)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Computes quaternion associated to an orthogonal
c               matrix by using Spurrier's algorithm (vector,scalar).

c      Inputs:
c         t(3,3) - Orthogonal matrix

c      Outputs:
c         q(4)   - Quaternion
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n(5),ii,j,k,l,ll,i
      real*8    xm,bigm,trace,dum, t(3,3),q(4)

      save

      data      n /1,2,3,1,2/

c     Computation of quaternion

      trace = t(1,1) + t(2,2) + t(3,3)
      ii    = 1
      xm    = t(1,1)

      if(t(2,2).gt.xm) then
        xm = t(2,2)
        ii = 2
      endif

      if(t(3,3).gt.xm) then
        xm = t(3,3)
        ii = 3
      endif

      bigm = max(trace,xm)

      if(bigm.eq.trace) then

        q(4) = 0.5d0*sqrt(1.d0 + trace)
        dum  = 0.25d0 / q(4)
        do i = 1,3
          j    = n(i+1)
          k    = n(i+2)
          q(i) = dum*(t(k,j) - t(j,k))
        end do ! i

      else

        q(ii) = sqrt(0.5d0*t(ii,ii) + 0.25d0*(1.d0 - trace))
        dum   = 0.25d0 / q(ii)
        j     = n(ii+1)
        k     = n(ii+2)
        q(4)  = dum*(t(k,j) - t(j,k))
        do ll = ii+1,ii+2
          l    = n(ll)
          q(l) = dum*(t(l,ii) + t(ii,l))
        end do ! ll

      endif

      end
