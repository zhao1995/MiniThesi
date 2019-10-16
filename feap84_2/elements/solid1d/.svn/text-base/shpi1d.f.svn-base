c$Id:$
      subroutine shpi1d(sg,xl,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Enhanced mode shaped functions

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'incshp.h'

      integer   ndm
      real*8    sg(2),xl(ndm,*)

c     Compute enhanced strain 'incompatible' shape functions

      shpi(1,1) = -sg(1)/(xl(1,2) - xl(1,1))
      shpi(2,1) =  0.0d0
      shpi(1,3) =  0.0d0
      shpi(2,3) =  1.d0

      end
