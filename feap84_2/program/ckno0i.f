c$Id:$
      logical function ckno0i(iv, nn )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check that an integer vector has a non-zero component

c      Inputs:
c         iv(*)  - Vector of integers
c         nn     - Length of vector

c      Outputs:
c         ckno0i - true of non-zero entries exist; else false
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,nn, iv(*)

      save

      do n = 1,nn
        if(iv(n).ne.0) then
          ckno0i = .true.
          return
         endif
      end do ! n

c     Set false to indicate vector is zero

      ckno0i = .false.

      end
