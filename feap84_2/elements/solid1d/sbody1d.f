c$Id:$
      subroutine sbody1d(d,xl, r,ndm,ndf ,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set lint for d(5), if 0 set to nel               27/03/2009
c       2. Use 'quard1d' and 'interp1d' for solution        01/04/2009
c       3. Add flag to call list                            03/03/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Body force computation

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'eldata.h'
      include   'elbody.h'
      include   'qudshp.h'

      integer    ndm, ndf, isw, l,i
      real*8     d(*), xl(ndm,*), r(ndf,*), bf(3)
      real*8     rr

c     Body force computations

      if(isw.eq.15) then
        bf(1) = bodyf(1)
      else
        call sbodyf(d, bf)
      endif

c     Set quadrature order

      call quadr1d(d)

c     Compute body loadings

      do l = 1,lint

        call interp1d(l, xl, ndm,nel,.false.)

        if(nint(d(16)).eq.3) then
          rr = 0.0d0
          do i = 1,nel
            rr = rr + shp1(2,i,l)*xl(1,i)
          end do ! i
          jac(l) = jac(l)*rr
        endif

c       Compute residual

        do i = 1,nel
          r(1,i)  = r(1,i) + shp1(2,i,l)*bf(1)*jac(l)
        end do ! j
      end do ! l

      end
