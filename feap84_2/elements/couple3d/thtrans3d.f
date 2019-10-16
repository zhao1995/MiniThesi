c$Id:$
      subroutine thtrans3d(d,xl,vl,s,r, nel,ndf,ndm,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    14/01/2010
c       1. Move 'j1' statement down near line 135           29/03/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute transient thermal effects for 3-d elements

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         vl(ndf,*) - Temperature rate for element
c         nel       - Number of element nodes
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Size of element arrays

c      Outputs:
c         s(nst,*)  - Consistent or interpolated mass
c         r(ndf,*)  - Element inertial force
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eltran.h'   ! ctan(2)

      logical   ttfl
      integer   nel,ndf,ndm,nst, ii,i1, jj,j1, l,lint
      real*8    xsj,dv, aj1,aj2,lfac,cfac
      real*8    d(*),xl(ndm,nel),vl(ndf,nel)
      real*8    s(nst,nst),r(ndf,nel)
      real*8    shp(4,64),sg(4,125),sv(5,16)

      save

c     Compute mass quadrature

      l = min(5,nint(d(5)))
      if(nel.eq.4) then
        if(nint(d(182)).gt.0) then
          call tint3dn(nel,lint,sv)
        else
          l =  2
          call tint3d (l,lint,sv)
        endif
        ttfl = .true.
      elseif(nel.eq.10) then
        l   =  14
        call tint3d (l,lint,sv)
        ttfl = .true.
      elseif(nel.eq.11) then
        if(nint(d(182)).gt.0) then
          call tint3dn(nel,lint,sv)
        else
          l =  14
          call tint3d (l,lint,sv)
        endif
        ttfl = .true.
      elseif(nel.eq.14 .or. nel.eq.15) then
        l   =  16
        call tint3d (l,lint,sv)
        ttfl = .true.
      else
        ttfl = .false.
        if(nint(d(182)).gt.0) then
          call int3dn(nel,lint,sg)
        else
          if(l.eq.0) then
            if(nel.le.8) then
              l = 2
            elseif(nel.le.27) then
              l = 3
            else
              l = 4
            endif
          endif
          call int3d(l,lint,sg)
        endif
      endif

c     Set heat capacity interpolation: consistent (1) and lumped (0)

      cfac = d(7)
      lfac = 1.d0 - cfac

c     Initialize heat capacity

      do l = 1,lint

c       Compute shape functions

        if(ttfl) then
          call tetshp(sv(1,l),xl,ndm,nel,xsj,shp)
          dv = sv(5,l)*xsj*d(4)*d(64)
        else
          call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)
          dv = sg(4,l)*xsj*d(4)*d(64)
        endif

c       Compute heat capacity

        j1 = 4
        do jj = 1,nel

c         Compute db = rho*c*shape*dv

          aj1      = shp(4,jj)*dv
          aj2      = cfac*aj1
          aj1      = lfac*aj1
          s(j1,j1) = s(j1,j1) + aj1
          i1 = 4
          do ii = 1,nel
            s(i1,j1) = s(i1,j1) + shp(4,ii)*aj2
            i1       = i1 + ndf
          end do ! ii
          j1 = j1 + ndf
        end do ! jj

      end do ! l

c     Compute inertial effect

      i1 = 4
      do ii = 1,nel
        j1 = 4
        do jj = 1,nel
          r(4,ii)  = r(4,ii) - vl(4,jj)*s(i1,j1)
          s(i1,j1) = s(i1,j1)*ctan(2)
          j1       = j1 + ndf
        end do ! jj
        i1 = i1 + ndf
      end do ! ii

      end
