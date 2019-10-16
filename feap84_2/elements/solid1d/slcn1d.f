c$Id:$
      subroutine slcn1d(sig,eps,shp,xsj,p,s,se,lint,nel,nes)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension shp(2,8,*)  (4 -> 8)                   02/04/2009
c       2. Dimension shp(2,20,*)  (8 -> 20) (for qudshp)    29/03/2011
c       3. Add and project strains                          01/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Project element variables to nodes

c      Inputs:
c        sig(nes,*) - Stresses at quadrature points
c        eps(3,*)   - Strains  at quadrature points
c        shp(2,8,*) - Shape functions at quadrature points
c        xsj(*)     - Volume element at quadrature points
c        lint       - Number of quadrature points
c        nel        - Number nodes on element
c        nes        - Dimension of stress array

c      Outputs:
c        p(nen)   - Weights for 'lumped' projection
c        s(nen,*) - Integral of variables
c        se(nen)  - Error projectors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'sdata.h'
      include  'prstrs.h'
      include  'strnum.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   nel,nes, i,l,lint
      real*8    p(*),s(nen,*),se(*),xsj(*) ,sig(nes,*),eps(3,*)
      real*8    shp(2,20,*), xg

      save

c     Lumped and consistent projection routine

      do l = 1,lint

c       Compute lumped projection and assemble stress integrals

        do i = 1,nel

          xg   = shp(2,i,l)*xsj(l)
          p(i) = p(i) + xg

c         Stress projections

          s(i,1) = s(i,1) + sig(1,l)*xg
          s(i,2) = s(i,2) + sig(2,l)*xg
          s(i,3) = s(i,3) + sig(3,l)*xg

c         Strain projections

          s(i,4) = s(i,4) + eps(1,l)*xg
          s(i,5) = s(i,5) + eps(2,l)*xg
          s(i,6) = s(i,6) + eps(3,l)*xg

c         Error estimation projection

          se(i)  = se(i)  + erav*xg

        end do ! i
      end do ! l

      iste = 6

      end
