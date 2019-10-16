c$Id:$
      subroutine ctiend(ics,neps,dnope,nope,ipos)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Renumbering of surfaces for tied nodes

c     Inputs:
c       ics(*)     - List of original nodes for each surface.

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i,ii,n,neps,dnope,nope
      integer  ics(dnope,neps),ipos(*)

      save

c     Do current surface renumbering: ii = ipos(ii) for node map

      do n = 1,neps
        do i = 1,nope
          ii = ics(i,n)
          if(ii.gt.0 .and. ipos(ii).ne.ii) then
            ics(i,n) = ipos(ii)
          endif
        end do ! i
      end do ! n

      end
