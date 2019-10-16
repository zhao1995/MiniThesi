c$Id:$
      subroutine ster1d(d,xl,ul,th,st,ndf,ndm,nel,nen)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set quadrature for d(5).eq.0                     26/03/2009
c       2. Use 'quard1d' and 'interp1d' for solution        01/04/2009
c       3. Add flag to call list                            03/03/2010
c       4. Dimension dd(6,6,2)                              06/09/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'adapt1.h'
      include  'adapt2.h'
      include  'errind.h'
      include  'hdata.h'
      include  'qudshp.h'
      include  'comblk.h'

      integer   ndf,ndm,nel,nen, i,j,l
      real*8    st(nen,*),xl(ndm,*),sig(6), epsp(3)
      real*8    sigp(3),d(*),eps(9,3),ul(ndf,*),th(*)
      real*8    dd(6,6,2),dr(6,6), xx,ta

      save

c     Error computation routine

      vfem   = 0.0d0
      vproj  = 0.0d0
      verror = 0.0d0
      vener  = 0.0d0
      venere = 0.0d0
      heta   = 0.0d0

c     Set quadrature formula

      call quadr1d(d)

      do l = 1,lint

        call interp1d(l, xl, ndm,nel,.false.)

c       Compute stresses

        call stra1d(d,xl,ul,th,shp1(1,1,l),ndf,ndm,nel,xx,ta,eps)
        call estrsd(d,ta,eps,sig,dd,dr)

        do i = 1,3
          sigp(i)= 0.0d0
        end do ! i

        do i = 1,nel
          do j = 1,3
            sigp(j) = sigp(j) + shp1(2,i,l)*st(i,j)
          end do ! j
        end do ! i

c       Compute projected strains

        call invert(dd,3,6)

        do i = 1,3
c         epsp(i) = ta
          epsp(i) = 0.0d0
        end do ! i

        do i = 1,3
          do j = 1,3
            epsp(i) = epsp(i) + dd(i,j,1)*sigp(j)
          end do ! j
        end do ! i

c       Compute integral of stress squares for error indicator use

        heta = heta + jac(l)

        do i = 1,3
          vfem   = vfem   + sig(i)*sig(i)*jac(l)
          vproj  = vproj  + sigp(i)*sigp(i)*jac(l)
          verror = verror + ((sigp(i)-sig(i))**2)*jac(l)
          vener  = vener  + sig(i)*eps(i,1)*jac(l)
          venere = venere + (sigp(i)-sig(i))*(epsp(i)-eps(i,1))*jac(l)
        end do ! i
      end do ! l

c     Set error indicators

      arsq   = arsq  + heta
      efem   = efem  + vfem
      eproj  = eproj + vproj
      eerror = eerror+ verror
      eener  = eener + vener
      eenere = eenere+ venere

      areai  = heta
      heta  =  d(50)*sqrt(heta)

      end
