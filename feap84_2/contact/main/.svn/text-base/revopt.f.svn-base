c$Id:$
      subroutine revopt(nrev,nrrev,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Construct reverse map of nodes

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none
      integer    n, nr,numnp, nrev(*),nrrev(*)

c     Construct reverse map

      do n = 1,numnp
        nr = nrev(n)
        nrrev(nr) = n
      end do ! n

      end
