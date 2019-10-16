c$Id:$
      subroutine modl1u(d,ta,eps,hn,h1,nh,istrt, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension dd(2) for use in some models           11/06/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: 1-d Constitutive Equation Update at start of time for
c              midpoint algorithms

c     Input parameters
c          d(*)    -  up to ndd-nud-1 material parameters
c          eps     -  strain at t_n+1
c          hn(nh)  -  history terms at t_n
c          h1(nh)  -  history terms at t_n+alpha
c          nh      -  number of history terms
c          istrt   -  Start state: 0 = elastic; 1 = last solution
c          isw     -  Element control parameter

c     Output parameters
c          h1(nh)  -  history terms at t_n+1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdat1.h'
      include  'ddata.h'
      include  'pmod2d.h'

      integer   nh,istrt, umat,uprm,isw
      real*8    ta, sig,dd(2), d(*),eps,hn(*),h1(*)

      save

      if(theta(3).ne.1.0d0) then

c       Extract analysis type: 1=plane stress; 2=plane strain; 3=axi

        uprm  = ndd - nud
        umat  = int(d(uprm)) - 100

c       Program material models

        if(umat.lt.0) then

          plasfl = int(d(40)).eq.1
          viscfl = int(d(40)).eq.2

c         P l a s t i c i t y

          if(plasfl) then

c           Interpolate plastic strain and update internal variables

            h1(2) = hn(2) + (h1(2) - hn(2))/theta(3)
            call plas1d(d,ta,eps,h1,istrt, sig,dd)

c         V i s c o e l a s t i c i t y

          elseif(viscfl) then

c           call visc1d(d,sig,h1(1),h1(5),gfac)

c         E l a s t i c i t y

          else

            call elas1d(d,ta,eps,sig,dd)

          end if

c       U s e r    M o d e l    I n t e r f a c e

        else

          call umod1d(umat,eps,ta,d(1),d(uprm+1),hn(1),h1(1),nh,1,istrt,
     &                sig,dd(1), isw)

        end if

      end if

      end
