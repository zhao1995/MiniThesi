c$Id:$
      subroutine bm2res(nlay,d,bz,hn,h1,nh,def, forc,aa, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute resultant model from layers

c      Inputs:
c        nlay      - Number of z-levels
c        d(*)      - Material parameters
c        bz(2,*)   - Z-coordinate and B-width
c        hn(*)     - History terms at t_n
c        h1(*)     - History terms at t_n+1
c        nh        - Number of history terms per level
c        def(3,2)  - Axial, shear, and bending strains

c      Outputs:
c        forc(3,2) - Force resultants: Axial, shear, bending
c        aa(3,3,2) - Modulus array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bm2com.h'
      include  'bm2str.h'
      include  'elplot.h'

      logical   rayl
      integer   i,ii,istrt,n,nlay,nn,nh,isw
      real*8    ta,eps(2),sig(2),dd(2),kg,ks(2), bb,zz,df,dz
      real*8    d(*),bz(2,*),hn(*),h1(*),def(3,2),forc(3,2),aa(3,3,2)

      save

      data     ta /0.0d0/

c     Set flag for Rayleigh damping

      rayl = (d(77).ne.0.0d0  .or. d(78).ne.0.0d0)

c     Set shear parameters to elastic

      kg    = d(1)*d(37)*0.5d0/(1.d0+d(2))
      ks(1) = kg*def(2,1)
      ks(2) = kg*def(2,2)

c     Initialize arrays

      do i = 1,3
        forc(i,1) = 0.0d0
        forc(i,2) = 0.0d0
        aa(i,1,1) = 0.0d0
        aa(i,2,1) = 0.0d0
        aa(i,3,1) = 0.0d0
        aa(i,1,2) = 0.0d0
        aa(i,2,2) = 0.0d0
        aa(i,3,2) = 0.0d0
      end do ! i

c     Compute constitution using Gauss-Lobbato quadrature in layers

      nn     = 1
      ii     = 1
      zz     = bz(1,1)
      eps(1) = def(1,1) - zz*def(3,1)
      eps(2) = def(1,2) - zz*def(3,2)
      istrt  = nint(d(84))
      call modl1d(d,ta,eps,hn(nn),h1(nn),nh,ii,istrt, dd,sig, isw)
      siglr(1) = sig(1)
      epslr(1) = eps(1)
      tt(7)    = sig(1)
      tt(8)    = eps(1)

      do n = 1,nlay-1

c       Layer halfthickness

        dz        = 0.5d0*(bz(1,n+1) - bz(1,n))

c       First point for each layer

        bb        = bz(2,n)*dz*sl(2,1)

        df        = bb*sig(1)
        forc(1,1) = forc(1,1) + df
        forc(2,1) = forc(2,1) + ks(1)*bb
        forc(3,1) = forc(3,1) - df*zz

        df        = bb*dd(1)
        aa(1,1,1) = aa(1,1,1) + df
        aa(1,3,1) = aa(1,3,1) - df*zz
        aa(3,3,1) = aa(3,3,1) + df*zz*zz
        aa(2,2,1) = aa(2,2,1) + kg*bb

        if(rayl) then
          df        = bb*sig(2)
          forc(1,2) = forc(1,2) + df
          forc(2,2) = forc(2,2) + ks(2)*bb
          forc(3,2) = forc(3,2) - df*zz

          df        = bb*dd(2)
          aa(1,1,2) = aa(1,1,2) + df
          aa(1,3,2) = aa(1,3,2) - df*zz
          aa(3,3,2) = aa(3,3,2) + df*zz*zz
          aa(2,2,2) = aa(2,2,2) + kg*bb

        endif ! rayl

c       Add remaining points

        do i = 2,int(d(102))

          zz = 0.5d0*((1.d0 - sl(1,i))*bz(1,n)
     &              + (1.d0 + sl(1,i))*bz(1,n+1))
          bb = 0.5d0*((1.d0 - sl(1,i))*bz(2,n)
     &              + (1.d0 + sl(1,i))*bz(2,n+1))*dz*sl(2,i)
          eps(1) = def(1,1) - zz*def(3,1)
          eps(2) = def(1,2) - zz*def(3,2)

c         Increment counters

          nn  = nn + nh
          ii  = ii + 1

c         Get constitution: Stress (sig) and moduli (dd)

          istrt  = nint(d(84))
          call modl1d(d,ta,eps,hn(nn),h1(nn),nh,ii,istrt, dd,sig, isw)

          siglr(ii) = sig(1)
          epslr(ii) = eps(1)
          tt(2*ii+5)= sig(1)
          tt(2*ii+6)= eps(1)

          df          = bb*sig(1)
          forc(1,1)   = forc(1,1) + df
          forc(2,1)   = forc(2,1) + ks(1)*bb
          forc(3,1)   = forc(3,1) - df*zz

          df          = bb*dd(1)
          aa(1,1,1)   = aa(1,1,1) + df
          aa(1,3,1)   = aa(1,3,1) - df*zz
          aa(3,3,1)   = aa(3,3,1) + df*zz*zz
          aa(2,2,1)   = aa(2,2,1) + kg*bb

          if(rayl) then
            df          = bb*sig(2)
            forc(1,2)   = forc(1,2) + df
            forc(2,2)   = forc(2,2) + ks(2)*bb
            forc(3,2)   = forc(3,2) - df*zz

            df          = bb*dd(1)
            aa(1,1,2)   = aa(1,1,2) + df
            aa(1,3,2)   = aa(1,3,2) - df*zz
            aa(3,3,2)   = aa(3,3,2) + df*zz*zz
            aa(2,2,2)   = aa(2,2,2) + kg*bb
          endif ! rayl

        end do ! i
      end do ! n

      aa(3,1,1) = aa(1,3,1)
      aa(3,1,2) = aa(1,3,2)

      tt(1)   = forc(1,1)
      tt(2)   = def (1,1)
      tt(3)   = forc(2,1)
      tt(4)   = def (2,1)
      tt(5)   = forc(3,1)
      tt(6)   = def (3,1)

      end
