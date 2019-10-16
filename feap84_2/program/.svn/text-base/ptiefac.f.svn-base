c$Id:$
      subroutine ptiefac(x,ix, ir, r1,r2, is, tol)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Tie 4-node elements with common faces

c      Inputs:
c         x(ndm,*)   - Nodal coordinates of mesh
c         ix(nen1,*) - Element connection list
c         r1         - Region/material 1 value
c         r2         - Region/material 2 value
c                       is = 0 for material comparison
c                       is = 1 for region   comparison
c         tol        - Gap tolerance on merge

c      Outputs:
c         ir(*)      - Merge list
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'

      integer    r1,r2, is, m,n

      integer    ix(nen1,numel), ir(numnp)
      real*8     x(ndm,numnp), tol

      save

      do n = 1,numel
        if(ix(nen1-is,n).eq.r1) then
          do m = 1,numel
            if(ix(nen1-is,m).eq.r2) then
              call ptieint(ix(1,n),ix(1,m),x, ir, tol)
            endif
          end do ! m
        endif
      end do ! n

      end
