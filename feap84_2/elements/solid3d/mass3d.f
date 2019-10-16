c$Id:$
      subroutine mass3d(d,xl,s,p,ndf,ndm,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 14/15 node tet computation                   06/07/2007
c       2. Change 'ord' to 'nel' on 'tetshp' calls          05/11/2007
c       3. Add direct call to quadr3d and interp3d          11/11/2008
c       4. Remove 'nel' from call to 'quadr3d'              23/01/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute mass matrix for 3-d brick elements

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Size of element arrays

c      Outputs:
c         s(nst,*)  - Consistent or interpolated mass
c         p(nst)    - Diagonal (lumped) mass
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'pmod2d.h'
      include  'qudshp.h'

      integer   ndf,ndm,nst, i,ii,i1, jj,j1, l
      real*8    d(*),xl(ndm,*),s(nst,nst),p(nst)
      real*8    aj1,lfac,cfac

      save

c     Compute mass matrix

      call quadr3d(d,.false.)

c     Set mass interpolation factor between consistent(1) & lumped(0)

      cfac = d(7)
      lfac = 1.d0 - cfac

      do l = 1,lint

c       Compute shape functions

        call interp3d(l, xl, ndm,nel)
        jac(l) = jac(l)*d(4)

c       Compute db = rho*shape*dv

        j1 = 0
        do jj = 1,nel

c         Compute lumped mass

          aj1 = shp3(4,jj,l)*jac(l)
          do i = 1,3
            p(j1+i)      = p(j1+i)      + aj1
            s(j1+i,j1+i) = s(j1+i,j1+i) + aj1*lfac
          end do ! i

c         Compute consistent mass matrix

          aj1 = aj1*cfac
          i1  = 0
          do ii = 1,nel
            do i = 1,3
              s(i1+i,j1+i) = s(i1+i,j1+i) + shp3(4,ii,l)*aj1
            end do ! i
            i1 = i1 + ndf
          end do ! ii
          j1 = j1 + ndf
        end do ! jj
      end do ! l

      end
