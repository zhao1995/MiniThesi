c$Id:$
      subroutine updlnk(u,nneq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set boundary displacement values from links

c     Inputs:
c         u(n,4)  - Incremental displacements at linked nodes
c         nneq    - Number of components

c     Outputs:
c         u(n,1)  - Displacements at nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    nneq, n
      real*8     u(nneq,*)

      do n = 1,nneq
        u(n,1) = u(n,1) + u(n,4)
        u(n,4) = 0.0d0
      end do ! n

      end ! updlnk
