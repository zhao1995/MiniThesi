c$Id:$
      subroutine xnumb( jp, neq, num, nb )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute number of blocks needed for out-of -core
c               direct solution by profile solver

c      Inputs:
c         jp(*)     - pointers for column of profile form.
c         neq       - number of active equations.
c         num       - maximum number of terms in each block.

c      Outputs:
c         nb        - number of blocks au(al) take on disk.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'debugs.h'
      include  'iofile.h'

      integer   neq, num, nb, ns,ne,nn
      integer   jp(neq)

      save

c     Partition au/al to equal size blocks

      nb        = 1
      ns        = 2
      nn        = 1

      do ne = 2,neq
        if(jp(ne)-jp(nn).ge.num) then
          nb         = nb + 1
          ns         = ne
          nn         = ne - 1
        endif
      end do ! ne

      if(debug) then
        write(iow,2000) nb
        if(ior.lt.0) then
          write(*,2000) nb
        endif
      endif

2000  format(5x,'Number Blocks =',i5)

      end
