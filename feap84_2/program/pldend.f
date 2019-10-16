c$Id:$
          subroutine pldend(ldfor,lddis,ldflg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    06/03/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: End load table and save data to files

c     Input:
c        ldfor  - Flags for force inputs
c        lddis  - Flags for displacement inputs
c        ldflg  - Flag indicating load table is open

c     Output:
c        ldflg  - Flag indicating load table is closed
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    ldfor, lddis, ldflg

c     Save any existing forces and displacements in LDNOD/LDVAL

      if(ldfor) then
        call pldseta(mr(np(265)), hr(np(26)), 1, ndf,numnp)
        ldfor = .false.
      endif
      if(lddis) then
        call pldseta(mr(np(265)), hr(np(26)+nneq), 2, ndf,numnp)
        lddis = .false.
      endif

c     Set load flag to false

      ldflg = .false.

      end
