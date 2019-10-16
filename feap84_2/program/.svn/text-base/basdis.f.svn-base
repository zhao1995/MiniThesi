c$Id:$
      subroutine basdis(fpro,w,dt,mp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set multiple support base displacement values

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'prlod.h'
      include  'prld1.h'

      integer   mp,n, fpro(*)
      real*8    d,v ,dt, pr, w(mp,*)

      save

      if(dt.gt.0.0d0) then

        do n = 1,mp
          d = w(n,1)
          v = w(n,2)
          if(fpro(n).gt.0) then
            pr = prldv(fpro(n))
          else
            pr = prop
          endif
          w(n,1) = pr
          w(n,2) = 2.d0/dt*(w(n,1) - d) - v
          w(n,3) = (w(n,2) - v)/dt
        end do ! n

      end if

      end
