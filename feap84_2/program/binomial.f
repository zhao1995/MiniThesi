c$Id:$
      subroutine binomial(nn, fact, bi)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate binomial coefficients and factorials

c      Inputs:
c         nn        - Order terms

c      Outputs:
c         fact(*)   - Factorials
c         bi(*)     - Binomial coefficients
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    nn,no, n
      real*8     bi(*), fact(*)

      no      = nn - 1
      fact(1) = 1.d0
      do n = 2,nn
        fact(n) = dble(n)*fact(n-1)
      end do ! n

      bi(1 ) = 1.d0
      bi(nn) = 1.d0
      do n = 1,nn-2
        bi(n+1) = fact(no)/(fact(no-n)*fact(n))
      end do ! n

      end
