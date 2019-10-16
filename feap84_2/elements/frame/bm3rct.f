c$Id:$
      subroutine bm3rct(d,hn,h1,nh,strain, stress,mhook, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 3-d constitution                             04/11/2008
c       2. Increase quadrature modulo 10 to 100 and arrays  15/11/2008
c          sy and sz to 10.  Use Legendre form for n > 6
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
      include 'elpdat.h'
      include 'elplot.h'
      include 'tdata.h'

      logical  onedf
      integer  ic,ii,istrt,jj,jmax,nrct,nn,nh,isw, nr,nqy,nqz,ny,nz
      real*8   ta,eps(3),sig(3),dd(3,3),sarea
      real*8   jiner,iner1,iner2,inerc,yy,zz,df,yl,zl,yr,zr,da
      real*8   d(*),hn(*),h1(*),strain(*),stress(*),mhook(6,6)
      real*8   sy(2,11),sz(2,11),ww

      save

      data     ta /0.0d0/

c     Compute constitution using Gauss-Lobbato quadrature in layers

      onedf = nint(d(188)).eq.0
      sarea = 0.d0
      iner1 = 0.d0
      iner2 = 0.d0
      inerc = 0.d0
      nn    = 1
      ii    = 1
      jj    = 0
      nrct  = int(d(101))
      do nr = 1, nrct
        nqy = nint(d(101+5*nr))/100
        nqz = mod(nint(d(101+5*nr)),100)
        jmax  = 2*(nqy + nqz) - 4  ! Hack for one block
        if(nqy.le.6) then
          call int1dl(nqy,sy)
        else
          call int1d (nqy,sy)
        endif
        if(nqz.le.6) then
          call int1dl(nqz,sz)
        else
          call int1d (nqz,sz)
        endif
        yl = d(97 +5*nr)
        zl = d(98 +5*nr)
        yr = d(99 +5*nr)
        zr = d(100+5*nr)

        da = 0.25d0*(yr-yl)*(zr-zl)
        do nz = 1, nqz
          zz = 0.5d0*((1.d0 - sz(1,nz))*zl + (1.d0 + sz(1,nz))*zr)
          do ny = 1, nqy
            yy = 0.5d0*((1.d0 - sy(1,ny))*yl + (1.d0 + sy(1,ny))*yr)
            ww = sy(2,ny)*sz(2,nz)*da
            eps(1) = strain(3) + zz*strain(4) - yy*strain(5)
            istrt  = nint(d(84))

c           1-d constitution

            if(onedf) then
              call modl1d(d,ta,eps,hn(nn),h1(nn),nh,ii,istrt, dd,sig,
     &                    isw)

              if(isw.eq.20) then
                if(min(ny,nz).eq.1 .or.ny.eq.nqy .or. nz.eq.nqz) then
                  ic = nint(elplt(1))
                  if(ic.eq.2 .or. ic.eq.3) then
                    df = h1(nn+ic-2-nh)
                  else
                    df = sig(1)
                  endif
                  jj = jj + 1
                  call bm3pcn(yy,zz,df,jj,jmax,nqy,nqz)
                end if ! perimeter test
              else
                sarea = sarea + ww
                iner1 = iner1 + ww*zz*zz
                iner2 = iner2 + ww*yy*yy
                inerc = inerc + ww*yy*zz

                tt(3*ii-2)= sig(1)
                tt(3*ii-1)= eps(1)
                tt(3*ii  )= dd(1,1)

                df        = sig(1)*ww
                stress(3) = stress(3) + df
                stress(4) = stress(4) + df*zz
                stress(5) = stress(5) - df*yy

                df        = dd(1,1)*ww
                mhook(3,3) = mhook(3,3) + df
                mhook(3,4) = mhook(3,4) + df*zz
                mhook(3,5) = mhook(3,5) - df*yy
                mhook(4,3) = mhook(4,3) + df*zz
                mhook(4,4) = mhook(4,4) + df*zz*zz
                mhook(4,5) = mhook(4,5) - df*zz*yy
                mhook(5,3) = mhook(5,3) - df*yy
                mhook(5,4) = mhook(5,4) - df*zz*yy
                mhook(5,5) = mhook(5,5) + df*yy*yy
              endif ! isw check

c           3-d constitution

            else
              eps(2) = strain(1) - zz*strain(6)
              eps(3) = strain(2) + yy*strain(6)
              call modlbm(d,ta,eps,yy,zz,ww,hn(nn),h1(nn),nh,ii,istrt,
     &                    mhook, stress, isw)

            endif ! onedf check
            ii        = ii + 1
            nn        = nn + nh
          end do !nz
        end do !ny
      end do !nr
      jiner  = iner1 + iner2

c     Set cross sectional properties

      if(d(32).le.0.0d0) then
        d(32) = sarea
      endif
      if(d(33).le.0.0d0) then
        d(33) = iner1
        d(34) = iner2
        d(35) = inerc
        d(36) = jiner
      endif

c     Compute shear and torsional stiffness values for 1-d case

      if(nint(d(188)).eq.0) then
        sarea      =  sarea*d(27)
        mhook(1,1) =  d(37)*sarea
        mhook(2,2) =  d(38)*sarea
        mhook(6,6) =  jiner*d(27)
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
