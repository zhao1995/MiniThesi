c$Id:$
      subroutine autoslid(slid,ic,jj,nslid)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Auto surface descriptions in 2-d

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    slid(5,*),ic(4),nslid, j,jj

      do j = 1,4
        slid(j,nslid) = ic(j)
      end do ! j
      slid(5,nslid) = jj

      end
