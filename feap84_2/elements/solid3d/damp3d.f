c$Id:$
      subroutine damp3d(d,xl,s,ndf,ndm,nst)

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
c      Purpose: Compute damping matrix for 3-d brick elements

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Size of element arrays

c      Outputs:
c         s(nst,*)  - Damping matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'pmod2d.h'
      include  'qudshp.h'

      integer   ndf,ndm,nst, i,ii,i1, jj,j1, l
      real*8    d(*),xl(ndm,*),s(nst,nst)
      real*8    dv, aj1

      save

c     Compute damping matrix

      call quadr3d(d,.false.)

c     Loop over quadrature points

      do l = 1,lint

c       Compute shape functions

        call interp3d(l, xl, ndm,nel)
        dv = jac(l)*d(70)

c       Compute db = rho*shape*dv

        j1 = 0
        do jj = 1,nel

c         Compute consistent damping matrix

          aj1 = shp3(4,jj,l)*dv
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
