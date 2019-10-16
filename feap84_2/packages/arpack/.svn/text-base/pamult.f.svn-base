c$Id:$
      subroutine pamult(ittyp,u,v,fp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/06/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Mass interface to Lanczos/Arnoldi solvers

c      Inputs:
c         ittyp  - Matrix storage type
c         u(*)   - Current vector
c         fp(*)  - Pointers for matrix

c      Outputs:
c         v(*)   - Matrix * u result
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'comblk.h'
      include   'p_int.h'

      integer    ittyp
      real*8     u(*),v(*)

c     In-core profile multiplication from factors

      if(ittyp.eq.-3) then
        call primul(hr(fp(3)),hr(fp(2)),hr(fp(1)),u,v,mr(fp(4)),neq)
      endif

      end
