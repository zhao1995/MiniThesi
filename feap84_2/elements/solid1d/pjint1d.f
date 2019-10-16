c$Id:$
      subroutine pjint1d(d,ul,th,shp,dvol,epsv,sig,p,ndf,ndm,lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension shp(2,8,*)  (4 -> 8)                   02/04/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: J-integral and Matrial Force Calculations for linear
c              elements

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*)  - Current nodal solution parameters
c         th         - Nodal temperatures
c         shp(2,8,*) - Shape functions
c         dvol(*)    - Volume element
c         epsv(*)    - Volumetric strain measure
c         sig(*,*)   - Stresses
c         ndf        - Degree of freedoms/node
c         ndm        - Mesh coordinate dimension
c         nst        - Element array dimension
c         lint       - Number quadrature points

c      Outputs:
c         p(ndf,*)   - Material force
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'eldata.h'
      include   'hdata.h'
      include   'prstrs.h'
      include   'comblk.h'

      integer    ndf,ndm,type, i,j,k,l, lint
      real*8     d(*),ul(ndf,nen),th(*),p(ndf,*)
      real*8     shp(2,8,*),dvol(*),sig(9,*),epsv(*)
      real*8     enmom(3,3),gradu(3,3),eps(6),thi, thv,ta, Weng

      data       eps / 6*0.0d0 /
      save

c     Set nodal potential values: Can be specified or computed

      type = nint(d(16))
      ta   = 0.0d0
      do l = 1,lint

        call grad1d(ul,shp(1,1,l),epsv(l),ndf,ndm,nel,gradu,type)

c       Compute strains

        eps(1) = gradu(1,1)
        eps(3) = gradu(3,3)

c       Compute contributions to J-integrals

        Weng = 0.5d0*(sig(1,l)*eps(1) + sig(3,l)*eps(3))

        enmom(1,1) = sig(1,l)*gradu(1,1) - Weng
        enmom(2,2) =                     - Weng
        enmom(3,3) = sig(3,l)*gradu(3,3) - Weng

c       Accumulate integral

        do k = 1,nel
          do i = 1,ndm
            thi = shp(i,k,l)*dvol(l)
            thv = thi*th(k)
            do j = 1,ndm
              j_int(j) = j_int(j) + enmom(i,j)*thv
              p(j,k)   = p(j,k)   + enmom(i,j)*thi
            end do ! j
          end do ! i
        end do ! k
      end do ! l

      end
