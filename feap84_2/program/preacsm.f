c$Id:$
      subroutine preacsm(r,ndf,numnp,iqpl,qpl )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Reaction sum

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none
      integer    ndf,numnp, iqpl(3), n, n1,n2,ii
      real*8     r(ndf,numnp), qpl

      save

      ii  = iqpl(1)
      n1  = max(1,min(iqpl(2),iqpl(3)))
      n2  = max(iqpl(2),iqpl(3))
      if(n2.eq.0) then
        n1 = 1
        n2 = numnp
      endif
      qpl = 0.0d0
      do n = n1,n2
        qpl = qpl - r(ii,n)
      end do ! n

      end
