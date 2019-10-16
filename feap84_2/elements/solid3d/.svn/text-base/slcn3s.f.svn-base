c$Id:$
      subroutine slcn3s(xl,sigl,dt,st,ndm,nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'strnum.h'

      integer   ndm, nel, ii, jj, kk
      real*8    dt(*),st(nen,*),xl(ndm,*),rst(3,8),xez(3)
      real*8    shp(4,8),sigl(10,*), sqr3, shpj, detj

      save

      data      rst/-1.d0,-1.d0,-1.d0,   1.d0,-1.d0,-1.d0,
     &               1.d0, 1.d0,-1.d0,  -1.d0, 1.d0,-1.d0,
     &              -1.d0,-1.d0, 1.d0,   1.d0,-1.d0, 1.d0,
     &               1.d0, 1.d0, 1.d0,  -1.d0, 1.d0, 1.d0/

      sqr3 = 1.d0/sqrt(3.d0)

c     Compute jacobian parts

      do ii = 1,8

c       Compute shape function * jacobian

        do jj = 1,3
          xez(jj) = sqr3*rst(jj,ii)
        end do ! jj

        call bjac3d ( xez , xl, ndm, shp, detj )

c       Compute stresses at nodes from history terms

        do jj = 1,nel
          shpj = shp(4,jj)*detj
          dt(jj) = dt(jj)     + shpj
          do kk = 1,6
            st(jj,kk) = st(jj,kk) + shpj*sigl(kk,ii)
          end do ! kk
          if(sigl(10,ii).ne.0.0d0) then
            st(jj,10) = st(jj,10) + shpj*sigl(10,ii)
          endif
        end do ! jj
      end do ! ii

      iste = 10

      end
