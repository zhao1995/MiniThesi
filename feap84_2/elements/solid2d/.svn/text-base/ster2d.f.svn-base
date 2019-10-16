c$Id:$
      subroutine ster2d(ix,d,xl,ul,th,st,ndf,ndm,nel,nen)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add direct call to quadr2d and interp2d          11/11/2008
c          Remove shp from argument.
c       2. Remove 'nel' from call to 'quadr2d'              23/01/2009
c       3. Dimension dd(5,5,2) to prevent access error      11/06/2011
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

      integer   ndf,ndm,nel,nen, i,j,ii, ix(*)
      real*8    st(nen,*),xl(ndm,*),sig(6), epsp(4)
      real*8    sigp(4),d(*),eps(9,3),ul(ndf,*),th(*)
      real*8    dd(6,6,2),dr(6,6), xx,yy, ta

      save

c     Error computation routine

      vfem   = 0.0d0
      vproj  = 0.0d0
      verror = 0.0d0
      vener  = 0.0d0
      venere = 0.0d0
      heta   = 0.0d0

c     Set quadrature formula

      call quadr2d(d,.true.)

      do ii = 1,lint

        call interp2d(ii, xl,ix,ndm,nel, .false.)

c       Compute stresses

        call strn2d(d,xl,ul,th,shp2(1,1,ii),ndf,ndm,nel,xx,yy,ta,eps)
        call estrsd(d,ta,eps,sig,dd,dr)

        do i = 1,4
          sigp(i)= 0.0d0
        end do ! i

        do i = 1,nel
          do j = 1,4
            sigp(j) = sigp(j) + shp2(3,i,ii)*st(i,j)
          end do ! j
        end do ! i

c       Compute projected strains

        call invert(dd,4,6)

        do i = 1,3
c         epsp(i) = ta
          epsp(i) = 0.0d0
        end do ! i
        epsp(4) = 0.0d0

        do i = 1,4
          do j = 1,4
            epsp(i) = epsp(i) + dd(i,j,1)*sigp(j)
          end do ! j
        end do ! i

c       Compute integral of stress squares for error indicator use

        heta = heta + jac(ii)

        do i = 1,4
          vfem   = vfem   + sig(i)*sig(i)*jac(ii)
          vproj  = vproj  + sigp(i)*sigp(i)*jac(ii)
          verror = verror + ((sigp(i)-sig(i))**2)*jac(ii)
          vener  = vener  + sig(i)*eps(i,1)*jac(ii)
          venere = venere + (sigp(i)-sig(i))*(epsp(i)-eps(i,1))*jac(ii)
        end do ! i
      end do ! ii

c     Set error indicators

      arsq   = arsq  + heta
      efem   = efem  + vfem
      eproj  = eproj + vproj
      eerror = eerror+ verror
      eener  = eener + vener
      eenere = eenere+ venere

      areai  = heta

c     Check for triangles

      if(nel.eq.3 .or. nel.eq.6 .or. nel.eq.7) then
        heta = heta + heta
      endif

      heta  =  d(50)*sqrt(heta)

      end
