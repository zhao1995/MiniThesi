c$Id:$
      logical function pmetis(domains,numedg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Rename p_metis to pmetis                         31/05/2007
c       2. Convert to cover function for smetis for
c          METIS v5.x                                       06/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Partition mesh using METIS

c     Inputs:
c        domains - Number of domains to create

c     Outputs:
c        numedg  - Unused; Output for Windows version only
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#      include   "cdata.h"
#      include   "pointer.h"
#      include   "comblk.h"


      integer    domains, numedg, smetis, rval

      save

      rval = smetis(numnp,mr(np(252)),mr(np(253)),domains,mr(np(254)))

      if (rval .eq. 1) then
        pmetis = .true.
      else
        pmetis = .false.
      endif


      end
