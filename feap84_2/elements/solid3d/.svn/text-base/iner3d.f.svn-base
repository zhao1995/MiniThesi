c$Id:$
      subroutine iner3d(d,xl,vl,al,s,r, nel,ndf,ndm,nst, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    21/12/2007
c       1. Add 'vl' for Rayleigh mass damping               28/12/2007
c       2. Increase arrays to store 64 node brick           03/02/2009
c       3. Set quadrature for d(5).eq.0                     26/03/2009
c       4. Increase arrays to store 125 node brick          20/12/2010
c       5. Recode for nurbs/t-splines                       01/01/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute inertial effects for 3-d elements
c               Includes effects of Rayleigh mass damping.

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         vl(ndf,*) - Velocity for element
c         al(ndf,*) - Acceleration for element
c         ctan3     - Mass tangent factor
c         nel       - Number of element nodes
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Size of element arrays

c      Outputs:
c         s(nst,*)  - Consistent or interpolated mass
c         r(ndf,*)  - Element inertial force
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eltran.h'   ! ctan(3)
      include  'qudshp.h'   ! shp3(4,*,*), jac(*), lint

      integer   nel,ndf,ndm,nst, isw, i, ii,jj, l
      integer   ia(125)
      real*8    dv,dvm, aj1,aj2,aj3,lfac,cfac
      real*8    d(*),xl(ndm,nel),vl(ndf,nel),al(ndf,nel)
      real*8    s(nst,nst),r(ndf,nel)
      real*8    ac(3)

      save

c     Set assembly pointer

      ia(1) = 0
      do i = 2,nel
        ia(i) = ia(i-1) + ndf
      end do ! i

c     Compute mass quadrature

      call quadr3d(d,.false.)

c     Set mass interpolation factor between consistent(1) & lumped(0)

      cfac = d(7)
      lfac = 1.d0 - cfac
      dvm  = ctan(3) + d(77)*ctan(2)

c     Do quadrature

      do l = 1,lint

c       Compute shape functions

        call interp3d(l, xl, ndm,nel)
        dv = jac(l)*d(4)

c       Compute acceleration

        do i = 1,3
          ac(i) = 0.0d0
          do ii = 1,nel
            ac(i) = ac(i) + shp3(4,ii,l)*(al(i,ii) + d(77)*vl(i,jj))
          end do ! ii
          ac(i) = ac(i)*cfac
        end do ! i

c       Compute mass

        do jj = 1,nel

c         Compute db = rho*shape*dv

          aj1 = shp3(4,jj,l)*dv
          aj2 = aj1*lfac

c         Inertial residual

          do i = 1,3
            r(i,jj) = r(i,jj) - ac(i)*aj1
     &                        - (al(i,jj) + d(77)*vl(i,jj))*aj2
          end do ! i

c         Compute inertial tangent

          if(isw.eq.3) then
            aj1 = aj1*dvm
            aj2 = cfac*aj1
            aj1 = lfac*aj1
            do i = 1,3
              s(ia(jj)+i,ia(jj)+i) = s(ia(jj)+i,ia(jj)+i) + aj1
            end do ! j
            do ii = 1,nel
              aj3 = shp3(4,ii,l)*aj2
              do i = 1,3
                s(ia(ii)+i,ia(jj)+i) = s(ia(ii)+i,ia(jj)+i) + aj3
              end do ! i
            end do ! ii
          endif
        end do ! jj

      end do ! l

      end
