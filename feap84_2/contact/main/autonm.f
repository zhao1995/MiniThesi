c$Id:$
      subroutine autonm(ix,ib,ip,ma,nen,nen1,numnp,numel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Set pointer to store slidelines

c      Inputs:

c      Outputs:
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    ma,nen,nen1,numnp,numel, ii,j,n

      integer    ix(nen1,numel), ib(numnp),ip(numnp)

      do n = 1,numnp
        ip(n) = 0
      end do ! n

      do n = 1,numel
        if(ma.eq.0 .or. ma.eq.ix(nen1,n)) then
          do j = 1,nen
            ii = ix(j,n)
            if(ii.gt.0) then
              if(ib(ii).gt.0) then
                ip(ii) = ip(ii) + 1
              endif
            endif
          end do ! j
        endif
      end do ! n

      do n = 2,numnp
        ip(n) = ip(n) + ip(n-1)
      end do ! n

      end
