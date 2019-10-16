c$Id:$
      subroutine damp1d(d,xl,s,ndf,ndm,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set quadrature for d(5).eq.0                     26/03/2009
c       2. Add call to quadr1d and interp1d                 02/04/2009
c       3. Add flag to call list                            03/03/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute damping matrix for plane & axisymmetric problem

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Size of element arrays

c      Outputs:
c         s(nst,*)  - Consistent or interpolated mass
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'pmod2d.h'
      include  'qudshp.h'

      integer   ndf,ndm,nst, ii,i1, jj,j1, l
      real*8    dv, aj1,xx
      real*8    d(*),xl(ndm,*),s(nst,nst)

      save

c     Compute damping matrix

      call quadr1d(d)

      do l = 1,lint

c       Compute shape functions

        call interp1d(l, xl, ndm,nel,.false.)
        dv = abs(jac(l))*d(70)
        if(stype.eq.3) then
          xx = 0.0d0
          do jj = 1,nel
            xx = xx + shp1(2,jj,l)*xl(1,jj)
          end do ! jj
          dv = dv*xx
        end if

c       Compute db = c*shape*dv

        j1 = 1
        do jj = 1,nel

c         Compute damping matrix

          aj1      = shp1(2,jj,l)*dv
          i1  = 1
          do ii = 1,nel
            s(i1,j1) = s(i1,j1) + shp1(2,ii,l)*aj1
            i1       = i1 + ndf
          end do ! ii
          j1 = j1 + ndf
        end do ! jj
      end do ! l

      end
