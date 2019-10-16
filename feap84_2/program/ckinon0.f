c$Id:$
      logical function ckinon0(v, nn )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    22/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check that real vector has a non-zero component

c      Inputs:
c         v(*)   - Vector of integer numbers
c         nn     - Length of vector

c      Outputs:
c         cknon0 - true of non-zero entries exist; else false.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,nn, v(*)

      save

      do n = 1,nn
        if(v(n).ne.0) then
          ckinon0 = .true.
          return
         endif
      end do ! n

c     Set false to indicate vector is zero

      ckinon0 = .false.

      end
