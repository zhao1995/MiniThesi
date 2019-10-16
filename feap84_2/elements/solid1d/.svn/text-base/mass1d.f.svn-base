c$Id:$
      subroutine mass1d(d,xl,s,p,ndf,ndm,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set quadrature to nel for d(5).eq.0              26/03/2009
c       2. Add call to quadr1d and interp1d                 02/04/2009
c       3. Add flag to call list                            03/03/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute mass matrix for plane and axisymmetric problems

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

      integer   ndf,ndm,nst, ii,i1, jj,j1, l
      real*8    d(*),xl(ndm,*),s(nst,nst),p(nst)
      real*8    dv, aj1,xx,lfac,cfac

      save

c     Set quadrature order

      call quadr1d(d)

c     Set mass factors

      cfac = d(7)
      lfac = 1.d0 - cfac

      do l = 1,lint

c       Compute shape functions

        call interp1d(l, xl, ndm,nel,.false.)
        dv = abs(jac(l))*d(4)
        if(stype.eq.3) then
          xx = 0.0d0
          do jj = 1,nel
            xx = xx + shp1(2,jj,l)*xl(1,jj)
          end do ! jj
          dv = dv*xx
        end if

c       Compute db = rho*shape*dv

        j1 = 1
        do jj = 1,nel

c         Compute lumped mass matrices

          aj1      = shp1(2,jj,l)*dv
          p(j1)    = p(j1)    + aj1
          s(j1,j1) = s(j1,j1) + aj1*lfac

c         Compute consistent mass matrix

          aj1 = aj1*cfac
          i1  = 1
          do ii = 1,nel
            s(i1,j1) = s(i1,j1) + shp1(2,ii,l)*aj1
            i1       = i1 + ndf
          end do ! ii
          j1 = j1 + ndf
        end do ! jj
      end do ! l

      end
