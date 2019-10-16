c$Id:$
      subroutine pldtabl(ldtab,ldnum,ldval,ldtyp,lsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/01/2009
c       1. Increase ldtab to store spin number of displ.    09/03/2009
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Inputs:
c       ldnum        - Load set number
c       ldval        - Load value
c       ldtyp        - Load type: 1 = force; 2 = displacement
c       lsw          - Load entry

c     Outputs:
c       ldtab(4,2,*) - Load table

c     Definitions:

c       ldtab(1,1,ldnum) - Force start pointer
c       ldtab(2,1,ldnum) - Force length
c       ldtab(3,1,ldnum) - Force proportional load number
c       ldtab(4,1,ldnum) - Force (unused)

c       ldtab(1,2,ldnum) - Displ start pointer
c       ldtab(2,2,ldnum) - Displ length
c       ldtab(3,2,ldnum) - Displ proportional load number
c       ldtab(4,2,ldnum) - Displ spin set number
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none
      integer    ldnum,ldval,ldtyp,lsw,ldtab(4,2,*)

      ldtab(lsw,ldtyp,ldnum) = ldval

      end
