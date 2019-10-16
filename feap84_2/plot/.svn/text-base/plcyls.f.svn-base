c$Id:$
      subroutine plcyls(x,sig, ncomp, comp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set 'dz' for radius = 0; modify ncomp order      19/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Compute a cylindrical component to plot

c     Inputs:
c        x(ndm,*)     - Nodal coordinates
c        sig(numnp,*) - Nodal coordinates
c        ncomp        - Number of component to plot
c                       1 = rr    4 = rz
c                       2 = zz    5 = zt
c                       3 = tt    6 = tr
c                      >6 = other (copy direct)

c     Output:
c        comp(*)      - Cylindrical component for ncomp
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'
      include   'pdata3.h'
      include   'pview.h'

      integer    ncomp,    n, i, j
      real*8     radius,   c2,  s2,  sc
      real*8     x(ndm,*), sig(numnp,*), comp(*), dz(2)

      save

c     In case 1-d problem

      dz(2) = 0.0d0
      j     = min(ndf,ncomp)

      do n = 1,numnp

c       Axial and large components unmodified

        if(ncomp.eq.3 .or. ncomp.gt.6) then

          comp(n) = sig(n,ncomp)

c       Other components transformed cylindrically

        else
          do i = 1,2
            dz(i) = x(i,n) - zview0(i)
          end do ! i
          radius = sqrt(dz(1)**2 + dz(2)**2)

          if(radius.gt.0.0d0) then
            do i = 1,2
              dz(i) = dz(i)/radius
            end do ! i
          else
            dz(1) = 1.0d0
            dz(2) = 0.0d0
          endif
          c2 = dz(1)*dz(1)
          s2 = dz(2)*dz(2)
          sc = dz(1)*dz(2)

c         Set radial and tangential components

          if    (ncomp.eq.1) then
            comp(n) =  c2*sig(n,1) + 2.d0*sc*sig(n,4) + s2*sig(n,2)
          elseif(ncomp.eq.2) then
            comp(n) =  s2*sig(n,1) - 2.d0*sc*sig(n,4) + c2*sig(n,2)
          elseif(ncomp.eq.4) then
            comp(n) = (c2 - s2)*sig(n,4) + sc*(sig(n,2) - sig(n,1))
          elseif(ncomp.eq.5) then
            comp(n) =  dz(2)*sig(n,5) - dz(1)*sig(n,6)
          elseif(ncomp.eq.6) then
            comp(n) =  dz(1)*sig(n,5) + dz(2)*sig(n,6)
          endif

        endif

      end do ! n

      end
