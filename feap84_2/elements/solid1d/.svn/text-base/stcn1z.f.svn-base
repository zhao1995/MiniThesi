c$Id:$
      subroutine stcn1z(xl,sig,shp,xsj,lint,ndm,nel,nen)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove debug mprint                              04/10/2008
c       2. Dimension shp(2,8,*)  (nen -> 8)                 02/04/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: ZZ Element Projection

c      Inputs:
c         xl(ndm,*)    - Element node coordinates
c         sig(nen,*)   - Element stresses at quadrature points
c         shp(2,8,*)   - Shape functions at quadrature points
c         xsj(*)       - Jacobian at quadrature point (time weight)
c         lint         - Number quadrature points
c         ndm          - Mesh spatial dimension
c         nel          - Number nodes on element
c         nen          - Dimension for stress and shape functions

c      Outputs:
c         ek(10,10)    - Contribution to projection matrix
c         est(30,10)   - Contribution for each projection component
c                        N.B. Returned in 'zzcom1.h'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'zzcom1.h'

      integer    ndm,nel,nen, i,j,l,lint, nps
      real*8     xsj(*),xl(ndm,*),shp(2,8,*),sig(nen,*),xx(4), cons

c     Compute zz-projection and assemble local stress integrals

      nps = 2
      if(nel.eq.3) then
        nps = 3
      elseif(nel.eq.4) then
        nps = 4
      end if

      xx(1) =  1.d0
      do l = 1,lint
        xx(2) = -xnodz(1)
        do i = 1,nel
          xx(2) = xx(2) + shp(2,i,l)*xl(1,i)
        end do ! i

c       Quadratic and cubic projections

        if(nps.ge.3) then
          xx(3) = xx(2)*xx(2)
        endif

c       Cubic projection

        if(nps.eq.4) then
          xx(4) = xx(3)*xx(2)
        endif

c       Accumulate matrix and projection contributions

        do i = 1,nps
          cons = xx(i)*xsj(l)
          do j = 1,nps
            ek(j,i)  = ek(j,i) + xx(j)*cons
          end do ! j
          do j = 1,3
            est(j,i) = est(j,i) + cons*sig(j,l)
          end do ! j
        end do ! i

      end do ! l

      end
