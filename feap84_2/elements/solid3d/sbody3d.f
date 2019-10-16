c$Id:$
      subroutine sbody3d(d,xl, r,ndm,ndf ,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 14/15 node tet computation                   06/07/2007
c       2. Change 'ord' to 'nel' on 'tetshp' calls          05/11/2007
c       3. Remove 'nel' from call to 'quadr3d'              23/01/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Body force computations

c      Inputs:
c        d(*)      - Material set parameters
c        xl(ndm,*) - Element nodal coordinates
c        ndm       - Mesh dimension
c        ndf       - Dofs/node
c        isw       - Option flag

c      Outputs:
c        r(ndf,*)  - Body force nodal vector
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'eldata.h'
      include   'elbody.h'
      include   'qudshp.h'

      integer    ndm, ndf, isw, l,i
      real*8     d(*), xl(ndm,*), r(ndf,*)
      real*8     body(3)

c     Body force computations

      if(isw.eq.15) then
        body(1) = bodyf(1)
        body(2) = bodyf(2)
        body(3) = bodyf(3)
      else
        call sbodyf(d, body)
      endif

c     Set quadrature order

      call quadr3d(d,.true.)

c     Compute body loadings

      do l = 1,lint

        call interp3d(l, xl, ndm,nel)

c       Compute residual

        do i = 1,nel
          r(1,i)  = r(1,i) + shp3(4,i,l)*body(1)*jac(l)
          r(2,i)  = r(2,i) + shp3(4,i,l)*body(2)*jac(l)
          r(3,i)  = r(3,i) + shp3(4,i,l)*body(3)*jac(l)
        end do ! i
      end do ! l

      end
