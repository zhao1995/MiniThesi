c$Id:$
      subroutine stcn2z(xl,sig,eps,shp,xsj,lint,ndm,nel,nen)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove debug print                               04/10/2008
c       2. Define element types more broadly                04/07/2011
c       3. Dimensioin shp to 64                             17/12/2013
c       4. Project to reduced quadrature points and add     01/01/2014
c          'eps' for projection of strains.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: ZZ (SPR) Element Projection

c      Inputs:
c         xl(ndm,*)    - Element node coordinates
c         sig(nen,*)   - Element stresses at quadrature points
c         eps(6,*)     - Element strains  at quadrature points
c         shp(3,nen,*) - Shape functions at quadrature points
c         xsj(*)       - Jacobian at quadrature point (time weight)
c         lint         - Number quadrature points
c         ndm          - Mesh spatial dimension
c         nel          - Number nodes on element
c         nen          - Dimension for stress and shape functions

c      Outputs:
c         ek(10,10)    - Contribution to projection matrix
c         est(30,10)   - Contribution for each projection component
c                        N.B. Returned in 'zzcom1.h'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'eldatp.h'
      include   'strnum.h'
      include   'zzcom1.h'

      integer    ndm,nel,nen, i,j,l,lint, linr, nps
      real*8     xsj(*),xl(ndm,*),shp(3,64,*),sig(nen,*),eps(6,*)
      real*8     sigm(6), epsm(6), jac(9), sg(2,9), nshp(16,9)
      real*8     xn(2,16), xx(10,9)
      real*8     cons

c     Compute zz-projection and assemble local stress integrals

      if(nel.le.4) then
        nps = 3                            ! Linear
      elseif(nel.gt.4 .and. nel.le.9) then
        nps = 6                            ! Quadratic
      elseif(nel.gt.9) then
        nps = 10                           ! Cubic
      endif

c     Compute coordinates at quadrature points

      do i = 1,lint
        xn(1:2,i) = 0.0d0
        do j = 1,nel
          xn(1:2,i) = xn(1:2,i) + shp(3,j,i)*xl(1:2,j)
        end do ! j
      end do ! i

c     Set reduced quadrature points

      call zquad2d(lint, linr, sg)

c     Loop over reduced quadrature points

      do l = 1,linr

c       Compute shape functions

        call nshp2d(sg(1,l), nshp(1,l), lint)

c       Project values to reduced points

        xx(1,l)   =  1.0d0
        sigm(1:6) =  0.0d0
        epsm(1:6) =  0.0d0
        jac(l)    =  0.0d0
        xx(2,l)   = -xnodz(1)
        xx(3,l)   = -xnodz(2)
        do i = 1,lint
          xx(2,l)   = xx(2,l)   + nshp(i,l)*xn(1,i)
          xx(3,l)   = xx(3,l)   + nshp(i,l)*xn(2,i)
          jac(l)    = jac(l)    + nshp(i,l)*xsj(i)
          sigm(1:6) = sigm(1:6) + nshp(i,l)*sig(1:6,i)
          epsm(1:6) = epsm(1:6) + nshp(i,l)*eps(1:6,i)
        end do ! i
        if(nps.gt.3) then      ! Cubic
          xx( 4,l) = xx(2,l)*xx(2,l)
          xx( 5,l) = xx(2,l)*xx(3,l)
          xx( 6,l) = xx(3,l)*xx(3,l)
        endif
        if(nps.eq.10) then      ! Cubic
          xx( 7,l) = xx(4,l)*xx(2,l)
          xx( 8,l) = xx(5,l)*xx(2,l)
          xx( 9,l) = xx(5,l)*xx(3,l)
          xx(10,l) = xx(6,l)*xx(3,l)
        endif

c       Accumulate matrix and stress/strain projections

        do i = 1,nps
          cons        = xx(i,l)*jac(l)
          ek(1:nps,i) = ek(1:nps,i) + cons*xx(1:nps,l)
          est( 1:6,i) = est( 1:6,i) + cons*sigm(1:6)
          est(7:12,i) = est(7:12,i) + cons*epsm(1:6)
        end do ! i

      end do ! l

      iste = 12       ! Number projected items

c     Do history plots if required

      if(hpltfl) then
        call hlcn2z(xx,jac,nshp,linr,lint,nps)
      end if

      end
