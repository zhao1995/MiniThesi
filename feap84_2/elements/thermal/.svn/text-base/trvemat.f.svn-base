c$Id:$
      subroutine trvemat(d, ta,gradt, hn,hn1,nh, flux,dd,rhoc, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    19/11/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: User Constitutive Model for Thermal MPI3 Exchanges

c     Input:
c          ta      -  Temperature
c          gradt(3)-  Thremal gradient
c          hn(*)   -  History terms at t_n
c          nh      -  Number history terms
c          isw     -  Solution option from element

c     Output:
c          d(*)    -  Material parameters
c          hn1(nh) -  History terms at point: t_n+1
c                     N.B. 1-d models use only sig(1)
c          dd(3,3) -  Current material tangent moduli
c          rhoc    -  Density times specific heat
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    nh,isw
      real*8     ta, rhoc
      real*8     d(*), gradt(3),hn(nh),hn1(nh), flux(*),dd(3,*)

c     Compute and output flux (q) and conductivities (moduli)

      save

c     DUMMY MODULE:  Multi-scale use only from 'openmpi'

      end
