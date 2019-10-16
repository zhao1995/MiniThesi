c$Id:$
      subroutine fbody3w(rhoa,omega,xl,ul, r,s, tanfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Compute rotational body residual and tangent for finite
c              displacement beam elements

c     Inputs:
c        rhoa       - Density times area
c        omega      - Angular velocity in radians/time
c        xl(ndm,*)  - Nodal coordinates
c        ul(ndf,*)  - Nodal displacements

c     Outputs:
c        r(ndf,nel) - Residual for body force
c        s(nst,nst) - Tangent  for body force
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'eldata.h'
      include   'sdata.h'

      logical    tanfl
      integer    lint, i,j,l, ii,i1, jj,j1

      real*8     rhoa, omega, xjac, btj
      real*8     xl(ndm,nel), ul(ndf,nel), r(ndf,nel), s(nst,nst)
      real*8     shp(2,3),sg(2,4), bf(3),bt(3,3), xx(3)

c     Set quadrature

      lint = min(4,nel+1)
      call int1d(lint, sg)

c     Do quadrature

      do l = 1,lint

c       Get shape functions

        call shp1d(sg(1,l),xl,shp,ndm,nel,xjac)
        xjac = xjac*sg(2,l)

c       Compute deformed coordinates

        do i = 1,ndm
          xx(i) = 0.0d0
          do j = 1,nel
            xx(i) = xx(i) + shp(2,j)*(xl(i,j) + ul(i,j))
          end do ! j
        end do ! i

c       Get body force and tangent vector

        call sbodyw(rhoa,omega,xx, bf,bt, tanfl)

c       Compute residual term

        do i = 1,nel
          btj = shp(2,i)*xjac
          do j = 1,ndm
            r(j,i) = r(j,i) + bf(j)*btj
          end do ! j
        end do ! i

c       Compute tangent

        if(tanfl) then
          j1 = 0
          do j = 1,nel
            do jj = 1,3
              do ii = 1,3
                btj = bt(ii,jj)*shp(2,j)*xjac
                i1 = 0
                do i = 1,nel
                  s(i1+ii,j1+jj) = s(i1+ii,j1+jj) + shp(2,i)*btj
                  i1 = i1 + ndf
                end do ! i
              end do ! ii
            end do ! jj
            j1 = j1 + ndf
          end do ! j
        endif

      end do ! l

      end
