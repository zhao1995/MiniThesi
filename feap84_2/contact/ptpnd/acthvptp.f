c$Id:$
      subroutine acthvptp (nset)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: ACTivate History VARiable 5

c      Purpose: Activation of contact variables requested by problem

c      Inputs :

c      Outputs:
c         nset    - # of history set required for contact pair
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_pair.h'

      logical   errck,active
      integer   nset

      save

      call cdebug0 ('  acthvptp',-1)

c     Activation of variables for surface match

      errck = active('nfacn',1)
      errck = active('iclsn',1)

c     Activation of variables for normal contact

      errck = active ('istgn',1)

      errck = active ('gapn' ,1)

c     Activation of variables for friction

      if (iffric.eq.1) then
        errck = active ('kts'  ,2)
      endif

c     Activation of Lagrange multiplier variables

      if(ifsolm.eq.2) then
        errck = active ('lagmn',2)
        errck = active ('lagmu',1)
      endif

c     Activation of augmentation variables

      if(ifaugm.gt.1) then
        errck = active ('augfn',1)
          if (iffric.eq.1) then
          errck = active ('augft',1)
        endif
      endif

c     Stop variable activation and define # of data sets

      errck = active('stop', nset)

      end
