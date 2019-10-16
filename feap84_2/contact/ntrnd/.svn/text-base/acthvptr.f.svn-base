c$Id:$
      subroutine acthvptr ()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:


c      Acronym: ACTivate History Variable contact element 1

c      Purpose: Activation of contact variables requested by problem

c      Inputs :

c      Outputs:
c         nset    - # of history set required for the contactpair
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_geom.h'
      include  'c_pair.h'

      logical   errck,active
      integer   nset

      save

      call cdebug0 ('  acthvptr',-1)

c     Activation of variables requested by combination of flags
c     Variables for normal contact

      errck = active ('istgn',1)
      errck = active ('gn',1)

c     Define Lagrange multiplier data

      if(ifsolm.eq.2) then
c       errck = active ('lagmn',2)
        errck = active ('lagmu',1)
      endif

c     Define Augmentation data

      if(ifaugm.gt.1) then
        errck = active ('augfn',1)
      endif

c     Define # of data sets required for the pair

      if(nope1.eq.1) then
        nset = neps1
      elseif(nope1.eq.4) then
        nset = neps1*4
      else
        nset = neps1+1
      endif

c     Stop variable activation defining # of data set

      errck = active('stop',nset)

      end
