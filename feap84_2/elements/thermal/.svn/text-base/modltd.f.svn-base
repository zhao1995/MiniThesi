c$Id:$
      subroutine modltd(d, ta,gradt, hn,hn1,nh, dd,flux,rhoc, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved
c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    19/11/2009
c       1. Use 'oelmt.h' to average density & specific heat 09/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Thermal material model driver

c     Input parameters
c          d(*)      -  up to ndd-nud-1 material parameters
c          ta        -  Temperature
c          gradt(3)  -  Temperature gradient

c          hn(*)     -  History terms at point (time = t_n)
c          nh        -  Number history tems at point

c          isw       -  FEAP switch value

c     Output parameters
c          hn1(*)    -  History terms at pint (time = t_n+1)
c          dd(3,3)   -  Thermal conductivity tensor
c          flux(3)   -  Thermal flux
c          rhoc      -  Density times specific heat
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'oelmt.h'
      include   'setups.h'

      integer    nh,isw, tmat
      real*8     ta, psi,cs,sn,c2,s2, rhoc
      real*8     d(*), gradt(3), hn(*),hn1(*), dd(3,3), flux(3)

c     Set thermal material type number

      tmat = nint(d(193))

c     Multiscale material model

      if(tmat.eq.1) then

        call trvemat(d, ta,gradt, hn,hn1,nh, flux,dd,rhoc, isw)

c     Fourier linear model

      elseif(tmat.eq.2) then

c       compute conductivity tensor

        psi = d(31)
        cs  = cos(psi)
        sn  = sin(psi)
        c2  = cs*cs
        s2  = sn*sn
        cs  = cs*sn

        dd(1,1) = c2*d(61) + s2*d(62)
        dd(1,2) = cs*(d(61) - d(62))
        dd(1,3) = 0.0d0

        dd(2,1) = dd(1,2)
        dd(2,2) = s2*d(61) + c2*d(62)
        dd(2,3) = 0.0d0

        dd(3,1) = 0.0d0
        dd(3,2) = 0.0d0
        dd(3,3) = d(63)

c       Set flux

        flux(1) = -dd(1,1)*gradt(1) - dd(1,2)*gradt(2)
        flux(2) = -dd(2,1)*gradt(1) - dd(2,2)*gradt(2)
        flux(3) = -dd(3,3)*gradt(3)

c       Integrate specific heat and density

        if(rank.gt.0) then
          v_avg = v_avg + el_vol
          v_rho = v_rho + el_rho
          v_c   = v_c   + el_c
        endif

c       Density * specific heat

        rhoc = d(4)*d(64)

      endif

      end
