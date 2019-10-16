c$Id:$
      subroutine rvematm(d,detf,f,gradt,ta, hn,hn1,nh, sig,flux,dd,kt,
     &                   isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    02/02/2010
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: RVE Thermo-mechanical exchange model

c     Input:
c          d(*)     -  Material parameters
c          detf     -  Determinant of deforamtion gradient
c          f(3,3)   -  Deformation gradient at pt (finite deformation)
c          gradt(3) -  Thermal gradient
c          ta       -  Temperature change
c          hn(nh)   -  History terms at point: t_n
c          nh       -  Number of history terms
c          isw      -  Solution option from element

c     Output:
c          hn1(nh)  -  History terms at point: t_n+1
c          sig(*)   -  Stresses at point.
c                      N.B. 1-d models use only sig(1)
c          flux(4)  -  Thermal flux and heat
c          dd(7,*)  -  Current material tangent moduli
c                      N.B. 1-d models use only dd(1,1) and dd(2,1)
c          kt(3,3)  -  Conductivity
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    nh,isw
      real*8     detf,ta
      real*8     d(*),f(3,3),gradt(3), hn(nh),hn1(nh), sig(*),flux(*)
      real*8     dd(7,*), kt(3,3)

c     Compute and output stress (sig) and (moduli)

      save

c     DUMMY module for RVE Thermo-mechanical

      end
