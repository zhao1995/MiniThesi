c$Id:$
      subroutine pjint2d(d,ul,th,shp,dvol,epsv,sig,p,ndf,ndm,lint,nes)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Fix dimension on shp( )                      SG  16/12/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: J-integral and Matrial Force Calculations for linear
c              elements

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*)  - Current nodal solution parameters
c         th         - Nodal temperatures
c         shp(3,*,*) - Shape functions
c         dvol(*)    - Volume element
c         epsv(*)    - Volumetric strain measure
c         sig(*,*)   - Stresses
c         ndf        - Degree of freedoms/node
c         ndm        - Mesh coordinate dimension
c         nst        - Element array dimension
c         lint       - Number quadrature points
c         nes        - Dimension for sig and shp arrays

c      Outputs:
c         p(ndf,*)   - Material force
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'eldata.h'
      include   'hdata.h'
      include   'prstrs.h'
      include   'comblk.h'

      integer    ndf,ndm,stype, i,j,k,l, lint, nes
      real*8     d(*),ul(ndf,nen),th(*),p(ndf,*)
      real*8     shp(3,64,*),dvol(*),sig(nes,*),epsv(*)
      real*8     enmom(3,3),gradu(3,3),eps(6),thi, thv,ta, Weng

      save

c     Set nodal potential values: Can be specified or computed

      stype = nint(d(16))
      ta    = 0.0d0
      do l = 1,lint

        call grad2d(ul,shp(1,1,l),epsv(l),ndf,ndm,nel,gradu,stype)

c       Compute strains

        eps(1) = gradu(1,1)
        eps(2) = gradu(2,2)
        eps(3) = gradu(3,3)
        eps(4) = gradu(1,2) + gradu(2,1)

c       Compute contributions to J-integrals

        Weng = 0.5d0*(sig(1,l)*eps(1) + sig(2,l)*eps(2)
     &              + sig(3,l)*eps(3) + sig(4,l)*eps(4))

        enmom(1,1) = sig(1,l)*gradu(1,1)
     &             + sig(4,l)*gradu(2,1) - Weng
        enmom(1,2) = sig(1,l)*gradu(1,2)
     &             + sig(4,l)*gradu(2,2)
        enmom(2,1) = sig(4,l)*gradu(1,1)
     &             + sig(2,l)*gradu(2,1)
        enmom(2,2) = sig(4,l)*gradu(1,2)
     &             + sig(2,l)*gradu(2,2) - Weng
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

      subroutine grad2d(ul,shp,epsv,ndf,ndm,nel,gradu,stype)

      implicit   none

      integer    ndf,ndm,nel, i,j,k, stype
      real*8     rr,epsv
      real*8     ul(ndf,*), shp(ndm+1,*), gradu(3,3)

c     Compute displacement gradient

      do i = 1,3
        do j = 1,3
          gradu(j,i) = 0.0d0
        end do ! j
      end do ! i
      do i = 1,ndm
        do j = 1,ndm
          do k = 1,nel
            gradu(j,i) = gradu(j,i) + ul(j,k)*shp(i,k)
          end do ! k
        end do ! j
      end do ! i

      if(stype.eq.3 .or. stype.eq.8) then
        rr         = 0.0d0
        gradu(3,3) = 0.0d0
        do k = 1,nel
          rr         = rr         + ul(1,k)*shp(3,k)
          gradu(3,3) = gradu(3,3) + ul(1,k)*shp(3,k)
        end do ! k

        if(rr.gt.0.0d0) then
          gradu(3,3) = gradu(3,3)/rr
        else
          gradu(3,3) = gradu(1,1)
        endif
      endif

c     Correct gradient for mixed model

      rr = (epsv - gradu(1,1) - gradu(2,2) - gradu(3,3))/3.d0

      gradu(1,1) = gradu(1,1) + rr
      gradu(2,2) = gradu(2,2) + rr
      gradu(3,3) = gradu(3,3) + rr

      end
