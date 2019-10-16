c$Id:$
      subroutine pedge4(nel,ix,xl,nen1,n)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Driver routine for pout4e

c      Inputs:
c         nel       - Number nodes on panel
c         ix(nen1,*)- List of nodes connected to elements
c         xl(3,*)   - Nodal coordinates for face
c         nen1      - Dimension of ix array
c         n         - Number of face to compute

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comblk.h'
      include  'pointer.h'

      integer   nel,nen1,n
      integer   ix(nen1,*)
      real*8    xl(3,*)

      save

      call pout4e(nel,ix,mr(np(64)),xl,hr(np(65)),mr(np(63)),nen1,3,n)

      end
