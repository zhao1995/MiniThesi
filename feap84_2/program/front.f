c$Id:$
      subroutine front(nfrnt,nd,n,numb,ntag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute front width for profile optimization.

c      Inputs:
c         n         - Number of node to add or remore
c         numb      - Number of nodes on front
c         ntag      - Indicator to add or remove node from front

c      Outputs:
c         nfrnt(*)  - Current nodes on front
c         nd(*)     - Indicator on active nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   l, m, n,numb,ntag, nfrnt(*),nd(*)

      save

c     Subroutine to update nodes on front

      if(ntag.eq.1) then

c       Add node to front if new node

        do m=1,numb
          if(n.eq.nfrnt(m)) return
        end do ! m
        numb = numb + 1
        nfrnt(numb) = n
        nd(n) = -1
      else

c       Remove node from front

        do m=1,numb
          if(n.eq.nfrnt(m)) go to 220
        end do ! m

  220   numb = numb - 1
        if(m.le.numb) then

          do l=m,numb
            nfrnt(l) = nfrnt(l+1)
          end do ! l
        endif
      endif

      end
