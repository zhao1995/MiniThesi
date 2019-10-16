c$Id:$
      subroutine baslod(pmass,w, mf,mp,i, bb)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set multiple support base mass forcel

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'modcon.h'

      integer  mf,mp,i,m
      real*8   bb
      real*8   pmass(mf,mp),w(mp,3)

      save

c     Compute Rayleigh damping factor

      do m = 1,mp
        bb = bb - pmass(i,m)*(w(m,3) + rayla0*w(m,2))
      end do ! m

      end
