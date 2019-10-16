c$Id:$
      subroutine prcylc(x,u, comp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c     2. Move dof > 2 into comp(i,*)                        10/10/2007
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Compute a cylindrical displacement to print

c     Inputs:
c        x(ndm,*)   - Nodal coordinates - reference
c        u(ndf,*)   - Nodal displacements

c     Output:
c        comp(ndf,*)    - Cylindrical components
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'
      include   'pconstant.h'
      include   'pview.h'

      integer    n, i
      real*8     x(ndm,*), u(ndf,*), comp(ndf,*), dz0(2)
      real*8     radius0,theta0

      save

      if(ndm.gt.1) then
        do n = 1,numnp

c         Radial and tangential components

          do i = 1,2
            dz0(i) = x(i,n) - zview0(i)
          end do ! i

c         Set radial and tangential components

          radius0   = sqrt(dz0(1)**2 + dz0(2)**2)
          theta0    = atan2( dz0(2), dz0(1))

          comp(1,n) = u(1,n)*cos(theta0) + u(2,n)*sin(theta0)
          comp(2,n) = u(2,n)*cos(theta0) - u(1,n)*sin(theta0)

c         Set axial and other values component

          do i = 3,ndf
            comp(i,n) = u(i,n)
          end do

c         Transform coordinates to cylindrical form

          if(cview) then
            theta0 = theta0*180.d0/pi
            if(theta0.lt.0.0d0) then
              theta0 = theta0 + 360.d0
            endif
            x(1,n) = radius0
            x(2,n) = theta0
            if(ndm.ge.3) then
              x(3,n) = x(3,n) - zview0(3)
            endif
          endif

        end do ! n
      endif

      end
