c$Id:$
      subroutine root (x,nn,dpn,pn1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006

c-----[--.----+-!--.----+----.----+------------------------------------]
c      Improve approximate root x; in addition we also obtain
c         dpn = derivative of p(n) at x
c         pn1 = value of p(n-1) at x
c-----[--.----+-!--.----+----.----+------------------------------------]
      implicit   none

      logical    notconv
      integer    nn, iter
      real*8     x,dpn,pn1, d,p,dp

      real*8     eps
      data       eps / 1.d-39 /

      iter = 0

      notconv = .true.
      do while(notconv .and. iter.lt.50)
        iter = iter + 1
        call recur (p,dp,pn1,x,nn)
        d = p/dp
        x = x - d
        if(abs(d).le.eps) then
          notconv = .false.
        endif
      end do ! while
      dpn = dp

      end
