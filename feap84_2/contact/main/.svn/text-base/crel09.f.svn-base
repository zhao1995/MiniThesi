c$Id:$
      subroutine crel09(nope,dnope, ix,es,ee, ixc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    09/01/2011
c       Store element type in dnope position                09/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set contact facet elements

c      Inputs:
c        nope         - Maximum number of nodes on ix/ixc facet
c        dnope        - Dimension for ixc array
c        ix(nen1,*)   - Element connection list from problem
c        es           - First element number to copy
c        ee           - Last element  number to copy

c      Outputs:
c        ixc(dnope,*) - Contact facet list
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'

      integer    nope,dnope, es,ee
      integer    ix(nen1,*), ixc(dnope,*)

      integer    i,n,ce

      ce = 0
      do n = es,ee
        ce = ce + 1
        do i = 1,nope
          ixc(i,ce) = ix(i,n)
        end do ! i
        ixc(dnope,ce) = ix(nen+6,n)
      end do ! n

      end
