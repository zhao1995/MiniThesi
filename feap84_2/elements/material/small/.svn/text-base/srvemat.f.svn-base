c$Id:$
      subroutine srvemat(d,eps,ta,hn,hn1,nh, sig,dd,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    10/12/2007
c       1. Add allocation of 'RVEMA' and set to ma          13/04/2009
c       2. Add temperature to sends (was umatl2)            13/06/2009
c       3. Add 'd' to argument list                         10/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: User Constitutive Model for MPI2 Exchanges

c     Input:
c          d(*)    -  Material parameters
c          eps(6)  -  Strain at point (small deformation)
c          ta      -  Temperature change
c          hn(nh)  -  History terms at point: t_n
c          nh      -  Number of history terms
c          isw     -  Solution option from element

c     Output:

c          hn1(nh) -  History terms at point: t_n+1
c          sig(*)  -  Stresses at point.
c                     N.B. 1-d models use only sig(1)
c          dd(6,*) -  Current material tangent moduli
c                     N.B. 1-d models use only dd(1,1) and dd(2,1)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    nh,isw
      real*8     ta
      real*8     d(*),eps(6),hn(nh),hn1(nh), sig(*),dd(6,*)

c     Compute and output stress (sig) and (moduli)

      save

c     DUMMY MODULE:  Multi-scale use only from 'openmpi'

      end
