c$Id:$
      subroutine bm3res(d,hn,h1,nh,strain, stress,mhook, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:
c        d(*)       - Material parameters
c        hn(*)      - History terms at t_n
c        h1(*)      - History terms at t_n+1
c        nh         - Number of history terms per level
c        strain(6)  - Axial, shear, and bending strains

c      Outputs:
c        stress(6)  - Force resultants: Axial, shear, bending
c        mhook(6,6) - Modulus array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'counts.h'
      include 'eldata.h'
      include 'elplot.h'
      include 'tdata.h'

      integer  i,nh,isw
      real*8   d(*),hn(*),h1(*),strain(*),stress(*),mhook(6,6)

      save

c     Initialize arrays

      do i = 1,6
        stress(i)  = 0.0d0
        mhook(i,1) = 0.0d0
        mhook(i,2) = 0.0d0
        mhook(i,3) = 0.0d0
        mhook(i,4) = 0.0d0
        mhook(i,5) = 0.0d0
        mhook(i,6) = 0.0d0
      end do ! i

c     Compute cross section radius, thickness, sector area

      if    (nint(d(100)).eq.1) then
        call bm3tub(d,hn,h1,nh,strain,stress,mhook, isw)
      elseif(nint(d(100)).eq.2) then
        call bm3rct(d,hn,h1,nh,strain,stress,mhook, isw)
      elseif(nint(d(100)).gt.2) then
        call bm3sec(d,hn,h1,nh,strain,stress,mhook,nint(d(100))-2, isw)
      end if

      end
