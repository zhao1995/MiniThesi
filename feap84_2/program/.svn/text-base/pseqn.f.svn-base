c$Id:$
      subroutine pseqn(ip,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set sequential numbers to integer vector

c      Inputs:
c         numnp  - Number nodes in mesh

c      Outputs:
c         ip(*)  - Equation numbers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n, numnp, ip(*)

      save

      do n = 1,numnp
        ip(n) = n
      end do ! n

      end
