c$Id:$
      subroutine stohman (hic,npair,cp0)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: STOre History variables MANagment

c      Purpose: Store  pointers for History variables

c      Inputs :
c         npair   - # of current pair

c      Outputs:
c         hic     - History Correspondence table
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_comnd.h'
      include  'c_keyh.h'

      integer   hic((c_lp1+c_lp3),*),npair, kv
      real*8    cp0(nr0,n0c3:nc03,*)

      save

      call cdebug0 ('  stohman',-1)

c     Store ch1 and ch3 offsets, set legth of the set, # of set

      cp0(2,-1,npair) = isgp1
      cp0(3,-1,npair) = isgp3
      cp0(4,-1,npair) = lh1
      cp0(5,-1,npair) = lh3
      cp0(6,-1,npair) = nset

c     Save memory location vectors

      do kv = 1,c_lp1
        hic(kv,npair) = p1(kv)
      end do
      do kv = 1,c_lp3
        hic(c_lp1+kv,npair) = p3(kv)
      end do

c     Update ch1 and ch3 offsets and total length

      isgp1 = isgp1 + lh1*nset
      isgp3 = isgp3 + lh3*nset
      tlch1 = isgp1 - 1
      tlch3 = isgp3 - 1

      end
