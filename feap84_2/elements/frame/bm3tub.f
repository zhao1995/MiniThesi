c$Id:$
      subroutine bm3tub(d,hn,h1,nh,strain, stress,mhook, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'constant.h'                    14/11/2006
c     2. Add 3-d constitution                               04/11/2008
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

      include 'bm2com.h'
      include 'bm2str.h'
      include 'counts.h'
      include 'eldata.h'
      include 'elplot.h'
      include 'pconstant.h'
      include 'tdata.h'

      logical  onedf
      integer  i,ii,istrt,j,ntub,nlob,nn,nh,isw
      real*8   ta,eps(3),sig(2),dd(2),thet1,thet2,theta,dthet,sarea
      real*8   garea,iner,rad,sn,cn,thick,yy,zz,df, ww
      real*8   d(*),hn(*),h1(*),strain(*),stress(*),mhook(6,6)

      save

      data     ta /0.0d0/


c     Compute cross section radius, thickness, sector area

      onedf = nint(d(188)).eq.0
      ntub  = int(d(101))
      nlob  = int(d(102))
      rad   = d(103)
      thick = d(104)
      sarea = pi*rad*thick
      dthet = 360.d0/dble(ntub)
      call int1dl(nlob,sl)

c     Set cross sectional  properties

      if(d(33).le.0.0d0) then
        iner  = sarea*rad*rad
        d(33) = iner
        d(34) = iner
        d(36) = iner + iner
      endif

      garea = sarea*2.d0
      if(d(32).le.0.0d0) then
        d(32) = garea
      endif

c     Compute constitution using Gauss-Lobbato quadrature in layers

      sarea = sarea/dble(ntub)
      nn    = 1
      ii    = 1
      thet1 = 0.d0
      do j = 1,ntub

c       Add remaining points

        thet2 = thet1 + dthet
        do i = 1,nlob - 1
          theta = 0.5d0*((1.d0 - sl(1,i))*thet1
     &                 + (1.d0 + sl(1,i))*thet2)
          call pdegree(theta, sn,cn)
          yy    = rad*cn
          zz    = rad*sn
          if(i.eq.1) then
            ww = sl(2,1) + sl(2,nlob)
          else
            ww = sl(2,i)
          endif
          eps(1) = strain(3) + zz*strain(4) - yy*strain(5)
          istrt = nint(d(84))

          if(onedf) then
            call modl1d(d,ta,eps,hn(nn),h1(nn),nh,ii,istrt, dd,sig, isw)

            siglr(ii) = sig(1)
            epslr(ii) = eps(1)
            tt(3*ii-2)= sig(1)
            tt(3*ii-1)= eps(1)
            tt(3*ii  )= dd(1)
            nout      = ii

            df        = sig(1)*sarea*ww
            stress(3) = stress(3) + df
            stress(4) = stress(4) + df*zz
            stress(5) = stress(5) - df*yy

            df        = dd(1)*sarea*ww
            mhook(3,3) = mhook(3,3) + df
            mhook(3,4) = mhook(3,4) + df*zz
            mhook(3,5) = mhook(3,5) - df*yy
            mhook(4,3) = mhook(4,3) + df*zz
            mhook(4,4) = mhook(4,4) + df*zz*zz
            mhook(4,5) = mhook(4,5) - df*zz*yy
            mhook(5,3) = mhook(5,3) - df*yy
            mhook(5,4) = mhook(5,4) - df*zz*yy
            mhook(5,5) = mhook(5,5) + df*yy*yy

c         3-d constitution

          else
            eps(2) = strain(1) - zz*strain(6)
            eps(3) = strain(2) + yy*strain(6)
            call modlbm(d,ta,eps,yy,zz,ww,hn(nn),h1(nn),nh,ii,istrt,
     &                  mhook, stress, isw)
          endif ! onedf
          ii        = ii + 1
          nn        = nn + nh

        end do !j
        thet1 = thet2
      end do !i

c     Compute shear and torsional stiffness values

      if(onedf) then
        garea      =  garea*d(27)
        mhook(1,1) =  d(37)*garea
        mhook(2,2) =  d(38)*garea
        mhook(6,6) =  garea*rad**2         ! 2 Pi r^3 t * G
        mhook(1,6) = -mhook(1,1)*d(95)
        mhook(6,1) =  mhook(1,6)
        mhook(2,6) =  mhook(2,2)*d(94)
        mhook(6,2) =  mhook(2,6)
        sarea      =  sarea/ntub

c       Compute shear and torsion

        stress(1)  =  mhook(1,1)*strain(1) + mhook(1,6)*strain(6)
        stress(2)  =  mhook(2,2)*strain(2) + mhook(2,6)*strain(6)
        stress(6)  =  mhook(6,1)*strain(1) + mhook(6,2)*strain(2)
     &             +  mhook(6,6)*strain(6)
      endif ! onedf

      end
