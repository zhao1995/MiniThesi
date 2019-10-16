c$Id:$
      subroutine hlcn2d(s,nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Project element history variables to nodes

c      Inputs:
c        nel          - Number nodes on element

c      Outputs:
c        s(nen,*)     - Integral of variables
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldatp.h'
      include  'qudshp.h'

      include  'pointer.h'
      include  'comblk.h'

      integer   ii, jj, l, nel
      real*8    s(nen,*), xj

      save

c     Initialize arrays

      do ii = 1,nel
        do jj = 1,plhmax
          s(ii,jj) = 0.0d0
        end do ! jj
      end do ! ii

c     Compute projections: int ( plhis * shp(i) * darea )

      do l = 1,lint
        do ii = 1,nel
          xj    = jac(l)*shp2(3,ii,l)
          do jj = 1,plhmax
            s(ii,jj) = s(ii,jj) + plhis(jj,l)*xj
          end do ! jj
        end do ! ii
      end do ! l

      end
