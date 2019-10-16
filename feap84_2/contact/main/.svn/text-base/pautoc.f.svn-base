c$Id:$
      subroutine pautoc(ma,ns)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Determine number of slidelines and their data.
c      Inputs:
c         ma      : Material number for search

c      Outputs:
c         ns      : Number of slideline surfaces found
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'iofile.h'
      include  'sdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   setvar, palloc
      integer   ma,ns

      save

      call cdebug0 ('  pautoc',-1)

c     Set the arrays to find patches

      setvar = palloc(136,'CTEM1', numnp*ndm , 2)
      setvar = palloc(137,'CTEM2', numnp     , 1)

c     Call routines to compute boundary patches

      call autonm(mr(np(33)),mr(np(78)),mr(np(137)),
     &            ma,nen,nen1,numnp,numel)

      call autocn(mr(np(33)),mr(np(78)),mr(np(137)),hr(np(206)),
     &            ma,ndm,nen1,ns)

      setvar = palloc(136,'CTEM1', 0,2)
      setvar = palloc(137,'CTEM2', 0,1)

      end
