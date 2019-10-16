c$Id:$
      subroutine slcn3du(sig,eps,shp,xsj, p,s, lint,nel,nes)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    02/02/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Project element variables to nodes for uniform stress
c               elements.

c      Inputs:
c        sig(10)      - Uniform stresses
c        shp(4,nes,*) - Shape functions at quadrature points
c        xsj(*)       - Volume element at quadrature points
c        lint         - Number of quadrature points
c        nel          - Number nodes on element
c        nes          - Dimension of shape function array

c      Outputs:
c        p(nen)       - Weights for 'lumped' projection
c        s(nen,*)     - Integral of variables
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'strnum.h'

      integer   ii, jj, l, lint, nel,nes
      real*8    p(*),s(nen,*), xsj(*),shp(4,nes,*),sig(*),eps(*), xj

      save

c     Initialize the arrays

      do ii = 1,nel
        p(ii)    = 0.0d0
        do jj = 1,12
          s(ii,jj) = 0.0d0
        end do ! jj
      end do ! ii

c     Compute projections: int ( sig * shp(i) * darea )

      do l = 1,lint
        do ii = 1,nel
          xj    = xsj(l)*shp(4,ii,l)
          p(ii) = p(ii)   + xj
          do jj = 1,6
            s(ii,jj  ) = s(ii,jj  ) + sig(jj)*xj
            s(ii,jj+6) = s(ii,jj+6) + eps(jj)*xj
          end do ! jj
        end do ! ii
      end do ! l

      iste = 12

      end
