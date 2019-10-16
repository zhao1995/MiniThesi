c$Id:$
      subroutine resdrb(mass,pi1,pin,thn,rcg,dt, prb)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute rigid body residual

c      Inputs:
c         mass     - Mass of rigid body
c         pi1      - Angular momentum at t_n+1
c         pin      - Angular momentum at t_n
c         thn      - Time interpolation value for t_n+thn
c         rcg(*)   - Solution/rate for rigid body
c         dt       - Time increment

c      Outputs:
c         prb(*)   - Residual
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'sdata.h'

      integer   i
      real*8    mass, dt, cons, cona, thn

      real*8    p(6),pi1(3), pin(3), rcg(3,11), prb(*)

      save

      cons = thn/dt
      cona = mass/dt

      do i=1,3
        p(i)   = - cona*(rcg(i,5) - rcg(i,10))
        p(i+3) = - cons*(pi1(i) - pin(i))
      end do ! i

c     2-D Case

      if(ndm.eq.2) then
        do i=1,2
          prb(i) = p(i)
        end do ! i
        prb(3) = p(6)

c     3-D Case

      else
        do i=1,6
          prb(i) = p(i)
        end do ! i
      endif

      end
