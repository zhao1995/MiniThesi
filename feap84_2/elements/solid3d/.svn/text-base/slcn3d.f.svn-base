c$Id:$
      subroutine slcn3d(sig,eps, p,s, nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add strain projections                           01/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Project element variables to nodes

c      Inputs:
c        sig(10,*)    - Stresses at quadrature points
c        eps( 6,*)    - Strains  at quadrature points
c        lint         - Number of quadrature points
c        nel          - Number nodes on element

c      Outputs:
c        p(nen)       - Weights for 'lumped' projection
c        s(nen,*)     - Integral of variables
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldatp.h'
      include  'qudshp.h'
      include  'strnum.h'

      include  'pointer.h'
      include  'comblk.h'

      integer   ii, jj, l, nel
      real*8    p(*),s(nen,*), sig(10,*), eps(6,*), xj

      save

c     Initialize arrays

      do ii = 1,nel
        p(ii)    = 0.0d0
        do jj = 1,12
          s(ii,jj) = 0.0d0
        end do ! jj
      end do ! ii

c     Compute projections: int ( sig * shp(i) * darea )

      do l = 1,lint
        do ii = 1,nel
          xj    = jac(l)*shp3(4,ii,l)
          p(ii) = p(ii)   + xj
          do jj = 1,6
            s(ii,jj  ) = s(ii,jj  ) + sig(jj,l)*xj
            s(ii,jj+6) = s(ii,jj+6) + eps(jj,l)*xj
          end do ! jj
        end do ! ii
      end do ! l

      iste = 12

c     Do history plots if required

      if(hpltfl) then
        call hlcn3d(hr(np(304)),nel)
      end if

      end
