c$Id:$
      subroutine bm3sec(d,hn,h1,nh,strain, stress,mhook,stype, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'constant.h'                  14/11/2006
c       2. Add 3-d constitution                             04/11/2008
c       3. Pass sig(1) to bm3scn                            09/01/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute beam resultant constitution for shaped sections:
c              stype = 1: Wide flange
c              stype = 2: Channel
c              stype = 3: Angle
c              stype = 4: Solid circular

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
      include 'pconstant.h'
      include 'tdata.h'

      logical  onedf
      integer  ii,nn,nh,stype,isw, nqudr, istrt
      real*8   ta,eps(3),sig(2),dd(2),sarea,xy(3,17),sw(3,17)
      real*8   d(*),hn(*),h1(*),strain(*),stress(*),mhook(6,6)
      real*8   jiner,iner1,iner2,inerc,yy,zz,df
      real*8   hh,wt,wb,tu,tb,tw,a1,a2,a3,ay,az, ww

      save

      data     ta /0.0d0/

c     Wide flange section

      onedf = nint(d(188)).eq.0
      if(stype.eq.1) then

        hh = d(101)
        wt = d(102)
        wb = d(103)
        tu = d(104)
        tb = d(105)
        tw = d(106)
        a1 = (hh-tu-tb)*tw
        a2 = wt*tu
        a3 = wb*tb
        ay = a1*(hh-tu+tb)*0.5d0 + a2*(hh-0.5d0*tu) + a3*tb*0.5d0
        ay = ay/(a1+a2+a3)

        xy(1, 1) = -0.50d0*wb
        xy(2, 1) = -ay
        xy(3, 1) =  0.25d0*a3

        xy(1, 2) =  0.50d0*wb
        xy(2, 2) = -ay
        xy(3, 2) =  0.25d0*a3

        xy(1, 3) =  0.50d0*wb
        xy(2, 3) =  tb-ay
        xy(3, 3) =  0.25d0*a3

        xy(1, 4) =  0.50d0*tw
        xy(2, 4) =  tb-ay
        xy(3, 4) =  0.25d0*a1

        xy(1, 5) =  0.50d0*tw
        xy(2, 5) =  hh-tu-ay
        xy(3, 5) =  0.25d0*a1

        xy(1, 6) =  0.50d0*wt
        xy(2, 6) =  hh-tu-ay
        xy(3, 6) =  0.25d0*a2

        xy(1, 7) =  0.50d0*wt
        xy(2, 7) =  hh-ay
        xy(3, 7) =  0.25d0*a2

        xy(1, 8) = -0.50d0*wt
        xy(2, 8) =  hh-ay
        xy(3, 8) =  0.25d0*a2

        xy(1, 9) = -0.50d0*wt
        xy(2, 9) =  hh-tu-ay
        xy(3, 9) =  0.25d0*a2

        xy(1,10) = -0.50d0*tw
        xy(2,10) =  hh-tu-ay
        xy(3,10) =  0.25d0*a1

        xy(1,11) = -0.50d0*tw
        xy(2,11) =  tb-ay
        xy(3,11) =  0.25d0*a1

        xy(1,12) = -0.50d0*wb
        xy(2,12) =  tb-ay
        xy(3,12) =  0.25d0*a3

        nqudr    =  12

c     Channel

      elseif(stype.eq.2) then

        hh = d(101)
        wt = d(102)
        wb = d(103)
        tu = d(104)
        tb = d(105)
        tw = d(106)
        a1 = (hh-tu-tb)*tw
        a2 = wt*tu
        a3 = wb*tb
        ay = a1*(hh-tu+tb)*0.5d0 + a2*(hh-0.5d0*tu) + a3*tb*0.5d0
        ay = ay/(a1+a2+a3)
        az = (a1*tw + a2*wt + a3*wb)*0.5d0
        az = az/(a1+a2+a3)

        xy(1, 1) = -az
        xy(2, 1) = -ay
        xy(3, 1) =  0.25d0*a3

        xy(1, 2) =  wb - az
        xy(2, 2) = -ay
        xy(3, 2) =  0.25d0*a3

        xy(1, 3) =  wb - az
        xy(2, 3) =  tb - ay
        xy(3, 3) =  0.25d0*a3

        xy(1, 4) =  tw - az
        xy(2, 4) =  tb - ay
        xy(3, 4) =  0.25d0*a1

        xy(1, 5) =  tw - az
        xy(2, 5) =  hh - tu - ay
        xy(3, 5) =  0.25d0*a1

        xy(1, 6) =  wt - az
        xy(2, 6) =  hh - tu - ay
        xy(3, 6) =  0.25d0*a2

        xy(1, 7) =  wt - az
        xy(2, 7) =  hh - ay
        xy(3, 7) =  0.25d0*a2

        xy(1, 8) = -az
        xy(2, 8) =  hh - ay
        xy(3, 8) =  0.25d0*a2

        xy(1, 9) = -az
        xy(2, 9) =  hh - tu - ay
        xy(3, 9) =  0.25d0*(a1 + a2)

        xy(1,10) = -az
        xy(2,10) =  tb - ay
        xy(3,10) =  0.25d0*(a1 + a3)

        nqudr    =  10

c     Angle

      elseif(stype.eq.3) then

        hh = d(101)
        wb = d(102)
        tw = d(103)
        tb = d(104)
        a1 = (hh-tu-tb)*tw
        a3 = wb*tb
        ay = (a1*(hh-tb) + a3*tb)*0.5d0
        ay = ay/(a1+a3)
        az = (a1*tw + a3*wb)*0.5d0
        az = az/(a1+a3)

        xy(1, 1) = -az
        xy(2, 1) = -ay
        xy(3, 1) =  0.25d0*a3

        xy(1, 2) =  wb - az
        xy(2, 2) = -ay
        xy(3, 2) =  0.25d0*a3

        xy(1, 3) =  wb - az
        xy(2, 3) =  tb - ay
        xy(3, 3) =  0.25d0*a3

        xy(1, 4) =  tw - az
        xy(2, 4) =  tb - ay
        xy(3, 4) =  0.25d0*a1

        xy(1, 5) =  tw - az
        xy(2, 5) =  hh - ay
        xy(3, 5) =  0.25d0*a1

        xy(1, 6) = -az
        xy(2, 6) =  hh - ay
        xy(3, 6) =  0.25d0*a1

        xy(1, 7) = -az
        xy(2, 7) =  tb - ay
        xy(3, 7) =  0.25d0*(a1 + a3)

        nqudr    =  7

c     Solid circular

      elseif(stype.eq.4) then

        hh    = d(101)             ! radius
        ii    = nint(d(102))       ! quadrature order
        a3    = pi*hh*hh  ! area weight

        call int2dc(ii,nqudr,sw)

        if(isw.eq.20) then
          if(nqudr.eq.17) then
            nqudr = 8
          else
            nqudr = 4
          endif
        endif
        do ii = 1,nqudr
          xy(1,ii) = hh*sw(1,ii)
          xy(2,ii) = hh*sw(2,ii)
          xy(3,ii) = a3*sw(3,ii)
        end do ! ii

      else
        write(*,*) ' BM3SEC: Section =',stype,' NOT CODED'
        call plstop()
      endif

c     Compute constitution using Gauss-Lobbato quadrature in layers

      sarea = 0.d0
      jiner = 0.d0
      iner1 = 0.d0
      iner2 = 0.d0
      inerc = 0.d0
      nn    = 1
      do ii = 1, nqudr
        yy = xy(1,ii)
        zz = xy(2,ii)
        ww = xy(3,ii)

        istrt  = nint(d(84))
        eps(1) = strain(3) + zz*strain(4) - yy*strain(5)

c       1-d constitution

        if(onedf) then
          call modl1d(d,ta,eps,hn(nn),h1(nn),nh,ii,istrt, dd,sig, isw)

          if(isw.eq.20) then
            call bm3scn(yy,zz,sig(1),ii,nqudr)
          else
            sarea = sarea + ww
            iner1 = iner1 + ww*zz*zz
            iner2 = iner2 + ww*yy*yy
            inerc = inerc + ww*yy*zz

            tt(3*ii-2) = sig(1)
            tt(3*ii-1) = eps(1)
            tt(3*ii  ) = dd(1)

            df         = sig(1)*ww
            stress(3)  = stress(3) + df
            stress(4)  = stress(4) + df*zz
            stress(5)  = stress(5) - df*yy

            df         = dd(1)*ww
            mhook(3,3) = mhook(3,3) + df
            mhook(3,4) = mhook(3,4) + df*zz
            mhook(3,5) = mhook(3,5) - df*yy
            mhook(4,3) = mhook(4,3) + df*zz
            mhook(4,4) = mhook(4,4) + df*zz*zz
            mhook(4,5) = mhook(4,5) - df*zz*yy
            mhook(5,3) = mhook(5,3) - df*yy
            mhook(5,4) = mhook(5,4) - df*zz*yy
            mhook(5,5) = mhook(5,5) + df*yy*yy
          endif ! isw

c       3-d constitution

        else
          eps(2) = strain(1) - zz*strain(6)
          eps(3) = strain(2) + yy*strain(6)
          call modlbm(d,ta,eps,yy,zz,ww,hn(nn),h1(nn),nh,ii,istrt,
     &                mhook, stress, isw)

        endif ! onedf
        nn        = nn + nh
      end do !ii
      jiner = iner1 + iner2

c     Set cross sectional area and inertias

      if(d(32).le.0.0d0) then
        d(32) = sarea
      endif
      if(d(33).le.0.0d0) then
        d(33) = iner1
        d(34) = iner2
        d(35) = inerc
        d(36) = jiner
      endif

c     Compute shear and torsional stiffness values

      if(onedf) then
        sarea      =  sarea*d(27)
        mhook(1,1) =  d(37)*sarea
        mhook(2,2) =  d(38)*sarea
        mhook(6,6) =  jiner*d(27)*0.2d0
        mhook(1,6) = -mhook(1,1)*d(95)
        mhook(6,1) =  mhook(1,6)
        mhook(2,6) =  mhook(2,2)*d(94)
        mhook(6,2) =  mhook(2,6)

c       Compute shear and torsion

        stress(1)  =  mhook(1,1)*strain(1) + mhook(1,6)*strain(6)
        stress(2)  =  mhook(2,2)*strain(2) + mhook(2,6)*strain(6)
        stress(6)  =  mhook(6,1)*strain(1) + mhook(6,2)*strain(2)
     &             +  mhook(6,6)*strain(6)
      endif

      end
