c$Id:$
      subroutine iner2d(d,xl,al,ix,s,r,ndf,ndm,nst, stype)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    27/12/2007
c       1. Set quadrature for d(5).eq.0                     26/03/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute inertial effects for 2-d elements

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         al(ndf,*) - Nodal acceleration for element
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Dimension of element arrays
c         stype     - 1,2 = Plane; 3,8 = Axisymmetric/+torsion

c      Outputs:
c         s(nst,*)  - Consistent or interpolated mass
c         p(nst)    - Element inertial force
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'eltran.h'
      include  'qudshp.h'

      integer   ndf,ndm,nst,stype, i,ii,i1, jj,j1, l, ndof
      real*8    cfac,lfac, dv,rr, aj1,aj2
      integer   ix(nel)
      real*8    d(*),xl(ndm,*),al(ndf,nel),s(nst,nst),r(ndf,nel)
      real*8    cmass(16,16)

      save

c     Set quadrature points and weights

      if(nel.eq.3) then
        if(d(182).gt.0.0d0) then
          call tint2dn(nel,lint,el2)
        else
          l = -3
          call tint2d(l,lint,el2)
        endif
        quad = .false.
      elseif(nel.eq.6 .or. nel.eq.7 ) then
        if(d(182).gt.0.0d0) then
          call tint2dn(nel,lint,el2)
        else
          l = 7
          call tint2d(l,lint,el2)
        endif
        quad = .false.
      else
        if(nint(d(182)).gt.0) then
          call int2dn(nel,lint,sg2)
        else
          l = min(5,nint(d(5)))
          if(l.eq.0) then
            if(nel.le.4) then
              l = 2
            elseif(nel.le.9) then
              l = 3
            else
              l = 4
            endif
          endif
          call int2d(l,lint,sg2)
        endif
        quad = .true.
      endif

c     Set mass interpolation factor between consistent (1) and lump (0)

      cfac = d(7)
      lfac = 1.d0 - cfac

c     Initalize mass

      do jj = 1,nel
        do ii = 1,nel
          cmass(ii,jj) = 0.0d0
        end do ! ii
      end do ! jj

c     Axisymmetry with torsion
      if(stype.eq.8) then
        ndof = 3
      else
        ndof = 2
      endif

c     Compute shape functions and derivatives in reference configuration

      do l = 1,lint
        if(quad) then
          call shp2d(sg2(1,l),xl,shp2,dv,ndm,nel,ix,.false.)
          dv = dv*sg2(3,l)*d(4)
        elseif(bsplfl) then ! 6-node B-Spline triangle
          call trispl(el2(1,l),xl,ndm, dv,shp2)
          dv = dv*el2(4,l)*d(4)
        else
          call trishp(el2(1,l),xl,ndm,nel-4,dv,shp2)
          dv = dv*el2(4,l)*d(4)
        endif

c       Compute axisymmetric radius

        if(stype.gt.2) then !
          rr = 0.0d0
          do i = 1,nel
            rr = rr + xl(1,i)*shp2(3,i,1)
          end do ! i
          dv = dv*rr
        endif

c       Compute mass

        do jj = 1,nel
          aj1 = shp2(3,jj,1)*dv
          aj2 = aj1*cfac
          aj1 = aj1*lfac
          cmass(jj,jj) = cmass(jj,jj) + aj1
          do ii = 1,nel
            cmass(ii,jj) = cmass(ii,jj) + shp2(3,ii,1)*aj2
          end do ! ii
        end do ! jj

      end do ! l

c     Compute inertial effect

      do ii = 1,nel
        do jj = 1,nel
          do i = 1,ndof
            r(i,ii) = r(i,ii) - al(i,jj)*cmass(ii,jj)
          end do ! i
        end do ! jj
      end do ! ii

c     Expand mass into element array

      j1 = 0
      do jj = 1,nel
        i1 = 0
        do ii = 1,nel
          cmass(ii,jj) = cmass(ii,jj) * ctan(3)
          do i = 1,ndof
            s(i1+i,j1+i) = s(i1+i,j1+i) + cmass(ii,jj)
          end do ! i
          i1 = i1 + ndf
        end do ! ii
        j1 = j1 + ndf
      end do ! jj

      end
