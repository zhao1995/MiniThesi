c$Id:$
      subroutine geomnd(shp,shpr,jac, sig, s, ndm,ndf,nel,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      1. Add coding to compute geometric stiffness           14/08/2010
c      2. Add axisymmetric terms                              27/08/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute geometric stiffness

c      Inputs:
c         shp(*,*)  - Shape function derivatives
c         shpr(*)   - Axisymmetric shape functions
c         jac       - Jacobian transformation
c         sig(*)    - Stress
c         ndm       - Mesh spatial dimension
c         ndf       - DOF/node (max)
c         nel       - Nodes on elements
c         nst       - Stiffness dimension

c      Outputs:
c         s(nst,*)  - Element geometric stiffness matrix
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      integer    ndm,ndf,nel,nst, i,j,k, i1,j1
      real*8     shp(ndm+1,*),shpr(*), jac, sig(*), s(nst,*)
      real*8     sigshp(2), fact

c     Compute stress gradient term

      j1 = 0
      do j = 1,nel
        sigshp(1) = (sig(1)*shp(1,j) + sig(4)*shp(2,j)
     &            +  sig(3)*shpr(j))*jac
        sigshp(2) = (sig(4)*shp(1,j) + sig(2)*shp(2,j))*jac
        i1 = 0
        do i = 1,nel
          fact = (shp(1,i) + shpr(i))*sigshp(1)
     &         +  shp(2,i)*sigshp(2)
          do k = 1,ndm
            s(i1+k,j1+k) = s(i1+k,j1+k) - fact
          end do ! k
          i1 = i1 + ndf
        end do ! i
        j1 = j1 + ndf
      end do ! j

      end
