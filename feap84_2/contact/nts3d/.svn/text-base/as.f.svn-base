c$Id:$
      subroutine as (kset,x,u,ix1,surpoin,aslave, node)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Anna Haraldsson             February 1998            1.0

c      Acronym: Area of Slave node.

c      Purpose: Compute area for slave nodes.

c      Inputs:
c         kset    - Slave node number
c         x(*)    - Nodal coordinate array
c         u(*)    - Nodal displacement array
c         ix1(*)  - Slave facet node connection array
c         surpoin - Surface points

c      Outputs:
c         aslave  - Area for node kset
c         node    - Global slave node number
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'ndata.h'
      include  'sdata.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   kset,ix1(dnope1,*),surpoin(*)
      integer   anzkn,node,anzseg,nel2,inod2,nod2,seg,j
      real*8    x(ndm,*),u(ndf,*),xm(3,4),a(3),b(3),c(3),aa,aslave

      fp(1)  = np(191) + surpoin(nsurf1)
      anzkn  = mr(fp(1))              ! number nodes on surface nsurf1
      node   = mr(fp(1)+kset)         ! kset-th node on surface nsurf1
      fp(1)  = mr(fp(1)+kset+anzkn)
      anzseg = mr(np(191)+fp(1))      ! size of patch of node kset

c     Quadrilateral facet surface area contribution

      if(nope1.eq.4) then
        aslave = 0.d0

        do seg = 1,anzseg
          nel2 = mr(np(191)+fp(1)+seg)  ! seg-th seg of patch node kset
          do nod2 = 1,4
            inod2 = ix1(nod2,nel2)
            do j = 1,3
              xm(j,nod2) = x(j,inod2) + u(j,inod2)
            end do ! j
          end do ! nod2

c         => xm(1..4) = nodes of seg-th segment of patch

          do j = 1,3
            a(j) = xm(j,4) - xm(j,2)
            b(j) = xm(j,1) - xm(j,2)
            c(j) = xm(j,3) - xm(j,2)
          enddo

c         Accumulate area

          aa     = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
          aslave = aslave
     &           + (sqrt((b(1)*b(1)+b(2)*b(2)+b(3)*b(3))*aa
     &                 - (a(1)*b(1)+a(2)*b(2)+a(3)*b(3))**2)
     &           +  sqrt((c(1)*c(1)+c(2)*c(2)+c(3)*c(3))*aa
     &                 - (a(1)*c(1)+a(2)*c(2)+a(3)*c(3))**2))*0.125d0

        end do ! seg

c     All other facet shapes

      else
        aslave = 1.0d0
      endif

      end
