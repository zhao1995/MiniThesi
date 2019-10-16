c$Id:$
      subroutine nshp2d(sg,nshp, lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/01/2014
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute 2-d shape functions

c      Inputs:
c        sg(2)    - Parent coordinates
c        lint     - Number element quadrature points

c      Outputs:
c        shp(*)   - Shape functions for point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   lint, k
      real*8    sg(2),nshp(*)

      integer   xi1(16),xi2(16)
      real*8    s1m,s1p, s2m,s2p, sp1,sp2, s12, sq1,sq2
      real*8    s1s2,s1s9, s2s2,s2s9, n1(4),n2(4)

      save

      data       xi1/1,2,2,1,3,4,2,2,4,3,1,1,3,4,4,3/
      data       xi2/1,1,2,2,1,1,3,4,2,2,4,3,3,3,4,4/

c     Factors

      s1m = 0.5d0 - 0.5d0*sg(1)
      s1p = 0.5d0 + 0.5d0*sg(1)
      s2m = 0.5d0 - 0.5d0*sg(2)
      s2p = 0.5d0 + 0.5d0*sg(2)

c     2 x 2  Quadrature
      if(lint.eq.4) then

        nshp(1) = s1m*s2m
        nshp(2) = s1p*s2m
        nshp(3) = s1p*s2p
        nshp(4) = s1m*s2p

c     3 x 3  Quadrature
      elseif(lint.eq.9) then

        s12     =  sg(1)*sg(2)
        sp1     =  s1m*s1p*4.d0
        sp2     =  s2m*s2p*4.d0
        nshp(1) =  s12*s1m*s2m
        nshp(2) = -s12*s1p*s2m
        nshp(3) =  s12*s1p*s2p
        nshp(4) = -s12*s1m*s2p
        nshp(5) = -sp1*s2m*sg(2)
        nshp(6) =  s1p*sp2*sg(1)
        nshp(7) =  sp1*s2p*sg(2)
        nshp(8) = -s1m*sp2*sg(1)
        nshp(9) =  sp1*sp2

c     4 x 4  Quadrature
      elseif(lint.eq.16) then

        sq1   = sg(1)*sg(1)
        sq2   = sg(2)*sg(2)
        s1s9  = (9.d0*sq1 - 1.d0)*0.0625d0
        s2s9  = (9.d0*sq2 - 1.d0)*0.0625d0
        s1s2  = (1.d0     -  sq1)*1.6875d0
        s2s2  = (1.d0     -  sq2)*1.6875d0

        n1(1)  = (1.d0 - sg(1))*s1s9
        n1(2)  = (1.d0 + sg(1))*s1s9
        n1(3)  = (one3 - sg(1))*s1s2
        n1(4)  = (one3 + sg(1))*s1s2

        n2(1)  = (1.d0 - sg(2))*s2s9
        n2(2)  = (1.d0 + sg(2))*s2s9
        n2(3)  = (one3 - sg(2))*s2s2
        n2(4)  = (one3 + sg(2))*s2s2

        do k = 1,16
          nshp(k) = n1(xi1(k))*n2(xi2(k))
        end do ! k

      endif

      end
