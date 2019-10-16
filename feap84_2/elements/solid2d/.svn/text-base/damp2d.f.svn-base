c$Id:$
      subroutine damp2d(d,xl,ix,s,ndf,ndm,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add direct call to quadr2d and interp2d          11/11/2008
c       2. Remove 'nel' from call to 'quadr2d'              23/01/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute damping matrix for plane/axisymmetric problems

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         ix(*)     - Element nodal connections
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

      integer   ndf,ndm,nst, i,ii,i1, jj,j1, l, ncp, ix(*)
      real*8    d(*),xl(ndm,*),s(nst,nst)
      real*8    aj1,xx

      save

c     Number of mass components

      if(stype.eq.8) then
        ncp = 3
      else
        ncp = 2
      endif

c     Compute mass matrix

      call quadr2d(d,.true.)

      do l = 1,lint

c       Compute shape functions

        call interp2d(l, xl,ix,ndm,nel, .true.)

        if(stype.eq.3 .or. stype.eq.8) then
          xx = 0.0d0
          do jj = 1,nel
            xx = xx + shp2(3,jj,l)*xl(1,jj)
          end do ! jj
          jac(l) = jac(l)*xx
        end if

c       Compute db = c*shape*dv

        j1 = 0
        do jj = 1,nel

c         Compute damping matrix

          aj1 = shp2(3,jj,l)*jac(l)
          i1  = 0
          do ii = 1,nel
            do i = 1,ncp
              s(i1+i,j1+i) = s(i1+i,j1+i) + shp2(3,ii,l)*aj1
            end do ! i
            i1 = i1 + ndf
          end do ! ii
          j1 = j1 + ndf
        end do ! jj
      end do ! l

      end
