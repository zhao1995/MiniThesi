c$Id:$
      subroutine pnumna(ix,nen1,nen,numel, ip)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pbody.h'

      integer   nen1,nen, numel, i,nn
      integer   ix(nen1,*), ip(*)

      save

c     Tag active nodes

      do nn = 1,numel
        if((ix(nen1-1,nn).ge.nreg1 .and. ix(nen1-1,nn).le.nreg2) .and.
     &     (maplt.eq.0 .or. ix(nen1,nn).eq.maplt) ) then
          do i = 1,nen
            if(ix(i,nn).gt.0) ip(ix(i,nn)) = 1
          end do ! i
        endif
      end do ! nn

c     Tag contact nodes

      call contact (318)

      end
