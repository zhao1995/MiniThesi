c$Id:$
      function dotid(a,b,id,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: dot (scalar) product of two vectors using id array.

c      Inputs:
c         a(*)  - Vector 1
c         b(*)  - Vector 2 to be accessed using id array
c         id(*) - Equation pointer array
c         nn    - length of vectors

c      Outputs:
c         dotid - Scalar product
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i,nn, id(*)
      real*8   dotid, a(*),b(*)

      save

      dotid = 0.d0
      do i = 1,nn
        if(id(i).gt.0) dotid = dotid + a(i)*b(id(i))
      end do ! i

      end
