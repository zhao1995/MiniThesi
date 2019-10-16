c$Id:$
      subroutine acthvt2 (nset)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor             6 March 2003            1.0

c      Acronym: ACTivate History VARiable for Tied 2d problems

c      Purpose: Activation of contact variables requested by problem

c      Inputs :

c      Outputs:
c         nset    - # of history set required for contact pair
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_geom.h'
      include  'c_pair.h'
      include  'sdata.h'

      logical   errck,active
      integer   nset

      save

      call cdebug0 ('  acthvt2',0)

c     Activation of variables for surface match

      errck = active ('mastl',2)

c     Activation of Lagrange multiplier variables

      if(ifsolm.eq.2) then
        errck = active ('lagm' ,nope1*ndm)
      endif

      errck = active('stop', nset)

      end
