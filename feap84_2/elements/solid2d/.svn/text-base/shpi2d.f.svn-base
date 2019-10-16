c$Id:$
      subroutine shpi2d(sg,xsj,xl,ndm)

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
      real*8    xsj, sg(2),xl(ndm,*)

      save

c     Compute enhanced strain 'incompatible' shape functions

      shpi(1,1) = -sg(1)*(-xl(2,1) - xl(2,2) + xl(2,3) + xl(2,4))/xsj
      shpi(2,1) =  sg(1)*(-xl(1,1) - xl(1,2) + xl(1,3) + xl(1,4))/xsj
      shpi(3,1) =  0.0d0
      shpi(1,2) =  sg(2)*(-xl(2,1) + xl(2,2) + xl(2,3) - xl(2,4))/xsj
      shpi(2,2) = -sg(2)*(-xl(1,1) + xl(1,2) + xl(1,3) - xl(1,4))/xsj
      shpi(3,2) =  0.0d0
      shpi(1,3) =  0.0d0
      shpi(2,3) =  0.0d0
      shpi(3,3) =  sg(1)*sg(2)

      end
