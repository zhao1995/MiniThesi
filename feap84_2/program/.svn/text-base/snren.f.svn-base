c$Id:$
      subroutine snren (numnp,nnn,nren)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Optimize profile numbering using algorithm of S. Sloan

c      Inputs:
c        numnp   - Number of nodes
c        nnn(*)  - Original numbering

c      Outputs:
c        nren(*) - Reverse numbering
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    n,numnp,nnn(*),nren(*)

c     Reverse numbering order for profile

      do n = 1,numnp
        nren(nnn(n)) = n
      end do ! n

      end
