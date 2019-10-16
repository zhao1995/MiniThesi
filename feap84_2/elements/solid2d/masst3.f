c$Id:$
      subroutine masst3(stype,cfac,rho,xl,al,p,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    29/03/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plane and axisymmetric 3-node triangle inertia terms

c      Inputs:
c         stype     - Axisymmetry indicator
c         cfac      - Consistent mass factor
c         rho       - Density
c         xl(ndm,*) - Element nodal coordinates
c         al(ndf,*) - Element nodal accelerations

c      Outputs:
c         p(ndf,*)  - Element intertia residual
c         s(ndf,*)  - Element mass terms
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'eltran.h'
      include   'sdata.h'

      logical    axfl, axtor
      integer    stype, i,i1,j,j1,l,lint,nd2
      real*8     rho,cfac, rr,dv
      real*8     xl(ndm,3),al(ndf,3),p(ndf,3),s(nst,nst)
      real*8     el(4,3),shp(3,3),ac(3)

      axfl  = stype.eq.3 .or.stype.eq.8
      axtor = stype.eq.8
      nd2   = ndf + ndf
      l     = -3                         ! 3-pt formula
      call tint2d(l,lint,el)
      do l = 1,lint
        call trishp(el(1,l),xl,ndm,1,dv,shp)
        dv = rho*el(4,l)*dv*cfac

c       Axisymmetric effect

        if(axfl) then
          rr = el(1,l)*xl(1,1) + el(2,l)*xl(1,2) + el(3,l)*xl(1,3)
          dv = dv*rr
        endif

c       Acceleration

        ac(1) = (shp(3,1)*al(1,1)
     &         + shp(3,2)*al(1,2)
     &         + shp(3,3)*al(1,3))*dv
        ac(2) = (shp(3,1)*al(2,1)
     &         + shp(3,2)*al(2,2)
     &         + shp(3,3)*al(2,3))*dv


c       Residual

        do j = 1,3
          p(1,j) = p(1,j) - shp(3,j)*ac(1)
          p(2,j) = p(2,j) - shp(3,j)*ac(2)
        end do ! j

c       Mass matrix

        dv             = dv*ctan(3)
        s(1    ,1    ) = s(1    ,1    ) + shp(3,1)*shp(3,1)*dv
        s(1    ,1+ndf) = s(1    ,1+ndf) + shp(3,1)*shp(3,2)*dv
        s(1    ,1+nd2) = s(1    ,1+nd2) + shp(3,1)*shp(3,3)*dv

        s(1+ndf,1    ) = s(1+ndf,1    ) + shp(3,2)*shp(3,1)*dv
        s(1+ndf,1+ndf) = s(1+ndf,1+ndf) + shp(3,2)*shp(3,2)*dv
        s(1+ndf,1+nd2) = s(1+ndf,1+nd2) + shp(3,2)*shp(3,3)*dv

        s(1+nd2,1    ) = s(1+nd2,1    ) + shp(3,3)*shp(3,1)*dv
        s(1+nd2,1+ndf) = s(1+nd2,1+ndf) + shp(3,3)*shp(3,2)*dv
        s(1+nd2,1+nd2) = s(1+nd2,1+nd2) + shp(3,3)*shp(3,3)*dv

c       Axisymmetric with torsion terms

        if(axtor) then
          ac(3) = (shp(3,1)*al(3,1)
     &           + shp(3,2)*al(3,2)
     &           + shp(3,3)*al(3,3))*dv
          do j = 1,3
            p(3,j) = p(3,j) - shp(3,j)*ac(3)
          end do ! j
        endif

      end do ! l

c     Fill in other mass terms

      j1 = 1
      do j = 1,3
        i1 = 1
        do i = 1,3
         s(i1+1,j1+1) = s(i1+1,j1+1) + s(i1,j1)
         if(axtor) then
           s(i1+2,j1+2) = s(i1+2,j1+2) + s(i1,j1)
         endif
         i1 = i1 + ndf
        end do ! i
        j1 = j1 + ndf
      end do ! j

      end
