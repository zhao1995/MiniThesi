c$Id:$
      subroutine sbody2d(d,xl,ix, r,ndm,ndf ,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add direct call to quadr2d and interp2d          11/11/2008
c       2. Remove 'nel' from call to 'quadr2d'              23/01/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Body force computation

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'eldata.h'
      include   'elbody.h'
      include   'qudshp.h'

      integer    ndm, ndf, isw, l,i, ix(*)
      real*8     d(*), xl(ndm,*), r(ndf,*), body(3)
      real*8     rr

c     Body force computations

      if(isw.eq.15) then
        body(1) = bodyf(1)
        body(2) = bodyf(2)
      else
        call sbodyf(d, body)
      endif

c     Set quadrature order

      call quadr2d(d,.true.)

c     Compute body loadings

      do l = 1,lint

        call interp2d(l, xl,ix, ndm,nel, .false.)

        if(nint(d(16)).eq.3) then
          rr = 0.0d0
          do i = 1,nel
            rr = rr + shp2(3,i,l)*xl(1,i)
          end do ! i
          jac(l) = jac(l)*rr
        endif

c       Compute residual

        do i = 1,nel
          r(1,i)  = r(1,i) + shp2(3,i,l)*body(1)*jac(l)
          r(2,i)  = r(2,i) + shp2(3,i,l)*body(2)*jac(l)
        end do ! j
      end do ! l

      end
