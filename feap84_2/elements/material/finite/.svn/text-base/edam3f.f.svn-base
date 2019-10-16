c$Id:$
      subroutine edam3f(d,wengy,xin,damg,ddam,load)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Damage variable function

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   load
      real*8    wengy, xin, damg, ddam, dex, hdam, xitr, xx, d(*)

c     Compute damage parameters

      xitr   = sqrt(2.d0*wengy)
      if(xitr .gt. xin) then
        load = .true.
        xin  = xitr
      endif

c     Compute damage function for xi-t(n+1)

      xx     = xin/d(2)
      if(xx.lt. 1.d-08) then
        hdam = 1.0d0 - xx*(0.5d0 - xx*(1.d0/6.d0 - xx/24.d0))
        ddam = -.5d0 + xx*(1.d0/3.d0 - xx*(.125d0 - xx/30.d0))
      else
        dex  = exp(-xx)
        hdam = (1.d0 - dex)/xx
        ddam = (dex - hdam)/xx
      endif

c     Set damage and tangent parameters

      damg   = d(1) + (1.d0 - d(1))*hdam
      ddam   = (1.d0 - d(1))/d(2)*ddam

      end
