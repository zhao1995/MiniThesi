c$Id:$
      subroutine pscalb(b,da,is,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Scaling of RHS and SOLUTION vector

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    is,neq, n
      real*8     b(*),da(*)

      save

      do n = is,neq
        b(n) = b(n)*da(n)
      end do ! n

      end
