c$Id:$
      subroutine mass2d(d,xl,ul,ix,s,p,ndf,ndm,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1.  Change 'nel.eq.4' to 'nel.eq.3'                 11/12/2006
c       2. Raise quadrature order for 3-node triangle       29/03/2007
c       3. Add direct call to quadr2d and interp2d          11/11/2008
c       4. Remove 'nel' from call to 'quadr2d'              23/01/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute mass matrix for plane and axisymmetric problems

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         ul(ndm,*) - Nodal displacements for element
c         ix(*)     - Element nodal connections
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

      integer   ndf,ndm,nst, i,ii,i1, jj,j1, l,ncp, ix(*)
      real*8    d(*),xl(ndm,*),ul(ndf,*),s(nst,nst),p(nst)
      real*8    aj1,xx,uu,lfac,cfac
      real*8    rfac(3)

      save

c     Number of mass components

      if(stype.eq.8) then
        ncp = 3
      else
        ncp = 2
      endif

c     Set default radius type

      do i = 1,ncp
        rfac(i) = 1.d0
      end do ! i

c     Compute mass matrix

      call quadr2d(d,.false.)

      cfac = d(7)
      lfac = 1.d0 - cfac

      do l = 1,lint

c       Compute shape functions

        call interp2d(l, xl,ix, ndm,nel, .true.)
        jac(l) = jac(l)*d(4)

        if(stype.eq.3 .or. stype.eq.8) then
          xx = 0.0d0
          do jj = 1,nel
            xx = xx + shp2(3,jj,l)*xl(1,jj)
          end do ! jj
          jac(l) = jac(l)*xx
        end if

c       Check for finite deformation torsion case

        if(dtype.ne.0 .and. stype.eq.8) then
          uu = 0.0d0
          do jj = 1,nel
            uu = uu + shp2(3,jj,l)*ul(1,jj)
          end do ! jj
          rfac(3) = (xx +uu)**2
        endif

c       Compute db = rho*shape*dv

        j1 = 0
        do jj = 1,nel

c         Compute lumped mass matrices

          aj1 = shp2(3,jj,l)*jac(l)
          do i = 1,ncp
            p(j1+i)      = p(j1+i)      + rfac(i)*aj1
            s(j1+i,j1+i) = s(j1+i,j1+i) + rfac(i)*aj1*lfac
          end do ! i

c         Compute consistent mass matrix

          aj1 = aj1*cfac
          i1  = 0
          do ii = 1,nel
            do i = 1,ncp
              s(i1+i,j1+i) = s(i1+i,j1+i) + rfac(i)*shp2(3,ii,l)*aj1
            end do ! i
            i1 = i1 + ndf
          end do ! ii
          j1 = j1 + ndf
        end do ! jj
      end do ! l

      end
