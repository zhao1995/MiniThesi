c$Id:$
      subroutine pflush(ifile)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    25/05/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Flush file buggers

c     Input:
c        ifile   - Logical unit number

c     Output:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    ifile

c     Flush buffer on units

      if(ifile.gt.0) then
        call flush(ifile)
      else
        call flush(ifile)
      endif

      end
