c$Id:$
      subroutine solgau (a,x,b,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nn, kc,kr,kl
      real*8    a(nn,nn),x(nn),b(nn), pivot

c     Triangularization

      do kc = 1,nn-1
        do kr = kc+1,nn
          if (a(kc,kc).ne.0.d0) then
            pivot = -a(kr,kc)/a(kc,kc)
            do kl = kc,nn
              a(kr,kl) = a(kr,kl) + a(kc,kl)*pivot
            end do
            b(kr) = b(kr) + b(kc)*pivot
          endif
        end do
      end do

c     Backsubstitution

      do kc = nn,2,-1
        do kr = kc-1,1,-1
          if (a(kc,kc).ne.0.d0) then
            pivot    = -a(kr,kc)/a(kc,kc)
            a(kr,kc) = a(kr,kc) + a(kc,kc)*pivot
            b(kr)    = b(kr)    + b(kc)   *pivot
          endif
        end do
      end do

c     Diagonal inversion

      do kr=1,nn
        if (a(kr,kr).ne.0.d0) then
          a(kr,kr) = 1.d0/a(kr,kr)
        endif
      end do

c     Multiply for right vector

      do kr=1,nn
        x(kr) = b(kr)*a(kr,kr)
      end do

      end
