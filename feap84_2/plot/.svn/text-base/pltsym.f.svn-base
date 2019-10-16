c$Id:$
      subroutine pltsym(x,ndm,numnp,nsy)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Reflect coordinates for symmetry plots

c      Inputs:
c         x(ndm,*)  - Coordinates to reflect
c         ndm       - Dimension of x array
c         numnp     - Number of nodes in mesh
c         nsy       - Symmetry type number to control which coord
c                     components to reflect

c      Outputs:
c         x(ndm,*)  - Reflected coordinates
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdatay.h'
      include  'pdatxt.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   n, i, ndm,numnp,nsy
      real*8    x(ndm,numnp)

      save

c     Set coordinates for symmetry

      fp(1) = npty - 1
      do n = 1, numnp
        if(mr(fp(1)+n).ge.0) then
          do i = 1,ndm
            x(i,n) = (x(i,n) - xsyc(i))*xsym(i,nsy) + xsyc(i)
          end do ! i
        endif
      end do ! n

      end
