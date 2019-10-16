c$Id:$
      subroutine acthvar1 ()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: ACTivate History VARiable 1

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

      call cdebug0 ('  acthvar1',-1)

c     Activation of variables requested by combination of flags
c     Variables for normal contact

      errck = active ('masts',1)
      errck = active ('lnc',1)
      errck = active ('istgt',1)
      errck = active ('istgn',1)
      errck = active ('s21',1)
      errck = active ('c21',1)
      errck = active ('d21',1)
      errck = active ('csi',1)
      errck = active ('gn',1)
      errck = active ('gt',1)
      errck = active ('area',1)
      errck = active ('s21c',1)
      errck = active ('c21c',1)
      errck = active ('csic',1)

      errck = active ('fn',1)
      errck = active ('dgnfn',1)

      errck = active ('fuoi',1)
      errck = active ('fnmax',1)
      errck = active ('fnresid',1)
      errck = active ('fnaug',1)

c     Variables for friction

      if (iffric.eq.1) then
        errck = active ('dtd',1)
        errck = active ('td',1)

        errck = active ('ft',1)
        errck = active ('dtdft',1)
        errck = active ('dgnft',1)
        errck = active ('tde',1)
        errck = active ('tdp',1)
        errck = active ('istfr',1)

        errck = active ('d21oi',1)
        errck = active ('csioi',1)
        errck = active ('csibb',1)
      endif

c     Variables for augmentation

      if (ifaugm.eq.2) then
        errck = active ('augfn',1)
        if (iffric.eq.1) then
          errck = active ('augft',1)
        endif

      elseif (ifaugm.eq.3) then
        errck = active ('augfn',1)
        errck = active ('istgn',3)
        errck = active ('fn',3)
        errck = active ('gn',3)
        errck = active ('daugfn',1)

      elseif ((ifaugm.eq.4) .or. (ifaugm.eq.5)) then
        errck = active ('augfn',1)
        errck = active ('augfn1',1)
        errck = active ('gn1',1)
        errck = active ('augfn2',1)
        errck = active ('gn2',1)

        if (iffric.eq.1) then
          errck = active ('augft',1)
          errck = active ('augft1',1)
          errck = active ('augft2',1)
          errck = active ('augstif1',1)
          errck = active ('augstif2',1)
        endif
      endif

c     Define Lagrange multiplier data

      if(ifsolm.eq.2) then
        errck = active ('lagmn',2)
        errck = active ('lagmu',1)
      endif

c     Define # of data sets required for the pair

      if(nope1.eq.1) then
        nset = neps1
      else
        nset = neps1+1
      endif

c     Stop variable activation defining # of data set

      errck = active('stop',nset)

      end
