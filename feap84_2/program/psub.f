c$Id:$
      function psub(val)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Accumulate subtraction of value

c     Inputs:
c        val   - Value to subtract

c     Outputs:
c        psub  - Value after subtraction or zero
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    psub, val, xval

      save

      data      xval / 0.0d0 /

c     Look at parameter

      if(val.eq.0.0d0) then
        xval = 0.0d0
      else
        xval = xval - val
      endif

      psub = xval

      end
