c$Id:$
      subroutine pjint3d(ul,th,epsv,sig,p,ndf,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: J-integral/Matrial Force Calculations of linear element

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'eldata.h'
      include   'hdata.h'
      include   'prstrs.h'
      include   'qudshp.h'
      include   'comblk.h'

      integer    ndf,ndm, i,j,k,l
      real*8     ul(ndf,nen),th(*),p(ndf,*),sig(10,*),epsv(*)
      real*8     enmom(3,3),gradu(3,3),eps(6),thi, thv,ta, Weng

      save

c     Set nodal potential values: Can be specified or computed

      ta   = 0.0d0
      do l = 1,lint

        call grad3d(ul,shp3(1,1,l),epsv(l),ndf,nel,gradu)

c       Compute strains

        eps(1) = gradu(1,1)
        eps(2) = gradu(2,2)
        eps(3) = gradu(3,3)
        eps(4) = gradu(1,2) + gradu(2,1)
        eps(5) = gradu(2,3) + gradu(3,2)
        eps(6) = gradu(3,1) + gradu(1,3)

c       Compute contributions to J-integrals

        Weng = 0.5d0*(sig(1,l)*eps(1) + sig(2,l)*eps(2)
     &              + sig(3,l)*eps(3) + sig(4,l)*eps(4)
     &              + sig(5,l)*eps(5) + sig(6,l)*eps(6))

        enmom(1,1) = sig(1,l)*gradu(1,1)
     &             + sig(4,l)*gradu(2,1)
     &             + sig(6,l)*gradu(3,1) - Weng
        enmom(1,2) = sig(1,l)*gradu(1,2)
     &             + sig(4,l)*gradu(2,2)
     &             + sig(6,l)*gradu(3,2)
        enmom(1,3) = sig(1,l)*gradu(1,3)
     &             + sig(4,l)*gradu(2,3)
     &             + sig(6,l)*gradu(3,3)
        enmom(2,1) = sig(4,l)*gradu(1,1)
     &             + sig(2,l)*gradu(2,1)
     &             + sig(5,l)*gradu(3,1)
        enmom(2,2) = sig(4,l)*gradu(1,2)
     &             + sig(2,l)*gradu(2,2)
     &             + sig(5,l)*gradu(3,2) - Weng
        enmom(2,3) = sig(4,l)*gradu(1,3)
     &             + sig(2,l)*gradu(2,3)
     &             + sig(5,l)*gradu(3,3)
        enmom(3,1) = sig(6,l)*gradu(1,1)
     &             + sig(5,l)*gradu(2,1)
     &             + sig(3,l)*gradu(3,1)
        enmom(3,2) = sig(6,l)*gradu(1,2)
     &             + sig(5,l)*gradu(2,2)
     &             + sig(3,l)*gradu(3,2)
        enmom(3,3) = sig(6,l)*gradu(1,3)
     &             + sig(5,l)*gradu(2,3)
     &             + sig(3,l)*gradu(3,3) - Weng

c       Accumulate integral

        do k = 1,nel
          do i = 1,ndm
            thi = shp3(i,k,l)*jac(l)
            thv = thi*th(k)
            do j = 1,ndm
              j_int(j) = j_int(j) + enmom(i,j)*thv
              p(j,k)   = p(j,k)   + enmom(i,j)*thi
            end do ! j
          end do ! i
        end do ! k
      end do ! l

      end

      subroutine grad3d(ul,shp,epsv,ndf,nel,gradu)

      implicit   none

      integer    ndf,nel, i,j,k
      real*8     ul(ndf,*), shp(4,*), gradu(3,3), rr,epsv

c     Compute displacement gradient

      do i = 1,3
        do j = 1,3
          gradu(j,i) = 0.0d0
          do k = 1,nel
            gradu(j,i) = gradu(j,i) + ul(j,k)*shp(i,k)
          end do ! k
        end do ! j
      end do ! i

c     Correct gradient for mixed model

      rr = (epsv - gradu(1,1) - gradu(2,2) - gradu(3,3))/3.d0

      gradu(1,1) = gradu(1,1) + rr
      gradu(2,2) = gradu(2,2) + rr
      gradu(3,3) = gradu(3,3) + rr

      end
