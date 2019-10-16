c$Id:$
      subroutine pextnd()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Dummy external node routine

c     Outputs:
c        NORMV array for nodal normals
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    setval,palloc

      setval = palloc(206,'NORMV',numnp*3,2)

      end
