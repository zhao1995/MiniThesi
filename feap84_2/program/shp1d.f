c$Id:$
      subroutine shp1d(s,xl,shp,ndm,nel,xjac)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute shape functions, natural derivatives, and
c              jacobian for 3-D beam at natural coordinate s.
c              Linear (2 node) or quadratic (3 node) element.
c                 o----------o                 o-----o-----o
c                 1          2                 1     3     2

c     Inputs:
c       s         : Natural coordinate
c       xl(3,nel) : Nodal global coordinates
c       ndm       : Coordinate dimension of mesh
c       nel       : Number of nodes of element

c     Outputs:
c       shp(2,nel): Shape functions and derivatives at s;
c                     shp(1,1 to nel): derivatives of shape functions
c                     shp(2,1 to nel): shape functions
c       xjac      : Jacobian at s
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   ndm,nel,i,j,k
      real*8    s,xjac, xl(ndm,nel),shp(2,nel)
      real*8    sp(10), xi,dxi,denom

      save

c     Linear element

      if(nel.eq.2) then

c       Computation of length of element

        xjac = 0.0d0
        do i = 1,ndm
          xjac = xjac + (xl(i,2) - xl(i,1))**2
        end do ! i
        if(xjac.eq.0.d0) then
           write(iow,3000) xjac
           call plstop()
        else
          xjac = sqrt(xjac)*0.5d0
        endif

        shp(2,1) = (1.d0 - s)*0.5d0
        shp(2,2) = (1.d0 + s)*0.5d0

        shp(1,1) = -0.5d0/xjac
        shp(1,2) =  0.5d0/xjac

c     Quadratic element

      elseif(nel.eq.3) then

c       Shape function values

        shp(2,1) =  s*(s - 1.d0)*0.5d0
        shp(2,2) =  s*(s + 1.d0)*0.5d0
        shp(2,3) =  1.d0 - s*s

c       Shape function natural derivatives

        shp(1,1) =  s - 0.5d0
        shp(1,2) =  s + 0.5d0
        shp(1,3) = -s*2.d0

c       jacobian

        xjac = 0.0d0
        do i = 1,ndm
          xjac = xjac + (shp(1,1)*xl(i,1)
     &                +  shp(1,2)*xl(i,2)
     &                +  shp(1,3)*xl(i,3))**2
        end do ! i
        if(xjac.eq.0.d0) then
           write(iow,3000) xjac
           call plstop()
        else
          xjac = sqrt(xjac)
        endif

c       Shape function global derivatives

        do i = 1,nel
          shp(1,i) = shp(1,i)/xjac
        end do ! i

c     Cubic to tenth order

      elseif(nel.le.10) then

c       Set location of nodes on parent element

        sp(1) = -1.0d0
        sp(2) =  1.0d0
        dxi   =  2.0d0/dble(nel - 1)
        xi    =  -1.d0
        do i = 3,nel
          xi    = xi + dxi
          sp(i) = xi
        end do ! i

c       Compute Lagrange interpolation and derivative

        do i = 1,nel
          shp(1,i) = 0.0d0
          shp(2,i) = 1.0d0
          denom    = 1.0d0
          do j = 1,nel
            if(i.ne.j) then
              shp(2,i) = shp(2,i)*(s  - sp(j))
              denom    = denom*(sp(i) - sp(j))
              dxi      = 1.0d0
              do k = 1,nel
                if(j.ne.k .and. i.ne.k) then
                  dxi = dxi*(s - sp(k))
                endif
              end do ! k
              shp(1,i) = shp(1,i) + dxi
            endif
          end do ! j
          shp(1,i) = shp(1,i)/denom
          shp(2,i) = shp(2,i)/denom
        end do ! i

c       Convert local derivatives to global ones

        xjac = 0.0d0
        do j = 1,ndm
          dxi = 0.0d0
          do i = 1,nel
            dxi = dxi + shp(1,i)*xl(j,i)
          end do ! i
          xjac = xjac + dxi*dxi
        end do ! j
        if(xjac.eq.0.d0) then
           write(iow,3000) xjac
           call plstop()
        else
          xjac = sqrt(xjac)
        endif

        do i = 1,nel
          shp(1,i) = shp(1,i)/xjac
        end do ! i

      endif

c     Format

3000  format(/5x,'Error in SHP1D: |J| =',1p,1d12.5/
     &       /5x,'Element has zero jacobian')

      end
