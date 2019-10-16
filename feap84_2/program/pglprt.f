c$Id:$
      subroutine pglprt(ug)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    12/04/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output global solution values

c      Inputs:
c         ug(*)  - Current solution values

c      Outputs:
c                - Global parameters output to files
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'bdata.h'
      include   'iofile.h'
      include   'pglob1.h'

      integer    n
      real*8     ug(*)

      if(ior.lt.0) then
        write(*,2000) head,(n,ug(n),n=1,geqnum)
      endif
      write(iow,2000) head,(n,ug(n),n=1,geqnum)

c     Format

2000  format(20a4/'   G l o b a l   E q u a t i o n   V a l u e s'/
     &       9x,'Equation No.',5x,'Value'/(7x,i10,1p,1e18.7))

      end
