c$Id:$
      subroutine pdxmsh(ndtyp,x)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Compute size of mesh in each coordinate direction.

c     Inputs:
c       ndtyp(*) - Node type: >=0 exist; < 0 tied
c       x(ndm,*) - Nodal coordinates

c     Outputs:
c       dxmsh(*) - Mesh size in each 'ndm' direction (via pdata0.h)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'pdata0.h'
      include   'sdata.h'

      logical    ndfl
      integer    ndtyp(*), i,n
      real*8     x(ndm,*), xmin,xmax

      save

      do i = 1,ndm
        ndfl = .true.
        do n = 1,numnp
          if(ndtyp(n).ge.0) then
            if(ndfl) then
              xmin = x(i,n)
              xmax = x(i,n)
              ndfl = .false.
            else
              xmin = min(xmin,x(i,n))
              xmax = max(xmax,x(i,n))
            endif
          endif
        end do ! n
        dxmsh(i) = max(1.d-8,xmax - xmin)
      end do ! i

      end
