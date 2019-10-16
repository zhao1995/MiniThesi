c$Id:$
      subroutine c2rplot(pen2,irsurf,crs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use pi from 'pconstant.h'                          14/11/2006
c     2. Set transformation data to plot rigid surface      15/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot Rigid surface

c      Inputs:
c        pen2   - Pen color to plot
c        irsurf - Surface definition
c        crs(:) - Data for surface

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pconstant.h'
      include   'pdata0.h'
      include   'pdatay.h'
      include   'prld1.h'
      include   'pview.h'
      include   'sdata.h'

      integer    pen2,irsurf, ip, n,nsy, j, dn,d1,d2
      real*8     crs(*),xc(3,4),xx(2),xg(3), rr,uu,pr, th,dth
      real*8     tr(3,3),xr(3)

c     Extract parameters

      dn = max(1,min(ndm,abs(nint(crs(1)))))
      rr = crs(2)
      uu = crs(3)*cs
      ip = nint(crs(4))

c     Set transformation data

      tr(1,1) = crs( 5)
      tr(2,1) = crs( 6)
      tr(3,1) = crs( 7)
      tr(1,2) = crs( 8)
      tr(2,2) = crs( 9)
      tr(3,2) = crs(10)
      tr(1,3) = crs(11)
      tr(2,3) = crs(12)
      tr(3,3) = crs(13)
      xr(1)   = crs(14)
      xr(2)   = crs(15)
      xr(3)   = crs(16)

c     Set proportional load factor

      if(ip.gt.0) then
        pr = prldv(ip)
      else
        pr = 0.0d0
      endif
      rr = rr + uu*pr

c     Cylindrical or Spherical surface

      if(irsurf.le.2) then

c       Plot surface

        call pppcol(pen2,0)

        th  = 0.0d0
        dth = pi/36.d0
        xx(1) =  rr
        do j = 1,3
          xg(j) = xr(j) + tr(j,1)*xx(1)
        end do ! j
        call plotl(xg(1),xg(2),xg(3),3)
        do n = 1,72
          th    = th + dth
          xx(1) = rr*cos(th)
          xx(2) = rr*sin(th)
          do j = 1,3
            xg(j) = xr(j) + tr(j,1)*xx(1) + tr(j,2)*xx(2)
          end do ! j
          call plotl(xg(1),xg(2),xg(3),2)
        end do ! n

c     Cartesian plane surface

      elseif(irsurf.eq.3) then

        do j = 1,4
          xc(dn,j) = rr
        end do ! j
        d1       = mod(dn,ndm) + 1
        xc(d1,1) = vmin(d1)
        xc(d1,2) = vmax(d1)
        xc(d1,3) = vmax(d1)
        xc(d1,4) = vmin(d1)
        if(ndm.eq.2) then
          do j = 1,4
            xc(3,j) = 0.0d0
          end do ! j
        elseif(ndm.eq.3) then
          d2       = mod(d1,ndm) + 1
          xc(d2,1) = vmin(d2)
          xc(d2,2) = vmin(d2)
          xc(d2,3) = vmax(d2)
          xc(d2,4) = vmax(d2)
        endif

c       Plot surface

        call pppcol(pen2,0)

c       Reflect coordinates for symmetry plots

        do nsy = 1,nsym
          lsym = isym(nsy)
          call pltsym(xc,3,4,lsym)
          if(ndm.eq.2) then
            call plotl(xc(1,1),xc(2,1),xc(3,1), 3)
            call plotl(xc(1,2),xc(2,2),xc(3,2), 2)
          elseif(ndm.eq.3) then
            call plotl(xc(1,4),xc(2,4),xc(3,4), 3)
            do j = 1,4
              call plotl(xc(1,j),xc(2,j),xc(3,j), 2)
            end do ! j
          endif
          call pltsym(xc,3,4,lsym)
        end do ! nsy

c     Other surface types

      else

        write(*,*) ' *WARNING* Other surface type plots not coded'

      endif

      end
