c$Id:$
      subroutine modl1d(d,td,eps,hn,h1,nh,ii,istrt, dd,sig,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Pass eps(1),sig(1) to elas1d                     09/01/2012
c       2. Pass eps(1),sig(1) to plas1d and visc1d          01/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Small Deformation 1-d Constitutive Equation Driver

c     Input parameters
c          d(*)    -  Up to ndd-nud-1 material parameters
c          td      -  Temperature
c          eps(2)  -  Current strain at point, and strain rate
c          hn(nh)  -  History terms at t_n
c          h1(nh)  -  History terms at t_n+1
c          nh      -  Number of history terms
c          ii      -  Number of calls to routine from each element
c          istrt   -  Start state: 0 = elastic; 1 = last solution
c          isw     -  Element control parameter

c     Ouput parameters
c          dd(2)   -  Current material tangent modulus
c          sig(2)  -  Stress at point.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdat1.h'
      include  'elcount.h'
      include  'pmod2d.h'

      integer   i,ii,istrt,nh, umat,uprm, isw
      real*8    td,theta(2), d(*),eps(*),hn(*),h1(*),dd(*),sig(*)

      save

c     Extract analysis type: 1=plane stress; 2=plane strain; 3=axi

      uprm  = ndd - nud
      umat  = int(d(uprm)) - 100

c     Set initial stress and zero tangent modulus

      dd(1)  = 0.0d0
      if(nint(d(160)).eq.1) then
        sig(1) = d(161)
      elseif(nint(d(160)).eq.2) then
      else
        sig(1) = 0.0d0
      endif

c     Program material models

      if(umat.lt.0) then

        plasfl = int(d(40)).eq.1
        viscfl = int(d(40)).eq.2

c       Move hn to h1

        do i = 1,nh
          h1(i) = hn(i)
        end do ! i

c       P l a s t i c i t y

        if(plasfl) then

          call plas1d(d,td,eps(1),h1,istrt, sig(1),dd)
          if(h1(3).eq.0.0d0) then
            nomats(1,2) = nomats(1,2) + 1
          else
            nomats(2,2) = nomats(2,2) + 1
          endif

c       V i s c o e l a s t i c i t y

        elseif(viscfl) then

          call visc1d(d,eps(1),h1(1),h1(2), sig(1),dd)
          nomats(1,4) = nomats(1,4) + 1

c       E l a s t i c i t y

        else

          call elas1d(d,td,eps(1), sig(1),dd)
          nomats(1,1) = nomats(1,1) + 1

        end if

c     U s e r    M o d e l    I n t e r f a c e

      else

        theta(1) = eps(1)
        theta(2) = eps(2)
        call umodel(umat,eps,theta,td,d,d(uprm+1),hn,h1,nh,
     &              ii,istrt,sig,dd, isw)

      end if

c     Check for tension/compression only materials

      if((d(167).gt.0.0d0 .and. sig(1).lt.0.0d0) .or.
     &   (d(167).lt.0.0d0 .and. sig(1).gt.0.0d0)) then
        sig(1) = 0.0d0
        dd(1)  = 0.0d0
      endif

c     Rayleigh stress rate

      sig(2) = dd(2)*eps(2)

      end
