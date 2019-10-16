c$Id:$
      subroutine plcylc(x,u, ncomp, comp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Compute a cylindrical component to plot

c     Inputs:
c        x(ndm,*)   - Nodal coordinates - reference
c        u(ndf,*)   - Nodal displacements
c        ncomp      - Number of component to plot
c                     1 = radial
c                     2 = tangential
c                     3 = axial

c     Output:
c        comp(*)    - Cylindrical component for ncomp
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'
      include   'pview.h'

      integer    ncomp,    n, i, j
      real*8     x(ndm,*), u(ndf,*), comp(*), dz0(2),dzn(2)
      real*8     radius0,radiusn

      save

c     In case 1-d problem

      dz0(2) = 0.0d0
      dzn(2) = 0.0d0
      j      = min(ndf,ncomp)

      do n = 1,numnp

c       Radial and tangential components

        if(ncomp.lt.3) then
          do i = 1,2
            dz0(i) = x(i,n) - zview0(i)
            dzn(i) = dz0(i) + u(i,n)
          end do ! i
          radius0 = sqrt(dz0(1)**2 + dz0(2)**2)
          radiusn = sqrt(dzn(1)**2 + dzn(2)**2)

c         Set radial and tangential components

          if(ncomp.eq.1) then
            comp(n) = radiusn - radius0
          else
            if(dz0(1).lt.0.0d0 .or. dzn(1).lt.0.0d0) then
              comp(n) = atan2(-dzn(2),-dzn(1))
     &                - atan2(-dz0(2),-dz0(1))
            else
              comp(n) = atan2( dzn(2), dzn(1))
     &                - atan2( dz0(2), dz0(1))
            endif
          endif

c       Set axial component

        else
          comp(n) = u(j,n)
        endif

      end do ! n

      end
