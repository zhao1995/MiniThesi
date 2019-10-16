c$Id:$
      subroutine upang(ia,ang,ul,ndf,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Update solution for sloping 2-d boundary

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    ia(*), ndf,isw, j
      real*8     ang,sn,cs,tm,ul(ndf,*)

      if(ia(1).ne.0 .and. ia(2).ne.0) then
        call pdegree(ang, sn,cs)
        if(isw.eq.1) then
          tm          =   cs*ul(ia(1),1) + sn*ul(ia(2),1)
          ul(ia(2),1) =  -sn*ul(ia(1),1) + cs*ul(ia(2),1)
          ul(ia(1),1) =   tm
        else
          do j = 1,3
            tm          =  cs*ul(ia(1),j) - sn*ul(ia(2),j)
            ul(ia(2),j) =  sn*ul(ia(1),j) + cs*ul(ia(2),j)
            ul(ia(1),j) =  tm
          end do ! j
        endif
      endif

      end
