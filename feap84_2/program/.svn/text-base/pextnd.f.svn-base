c$Id:$
      subroutine pextnd()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute external nodes on mesh

c      Inputs:
c         None

c      Outputs:
c         External nodes stored in mr(np(78))
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'qudshp.h'
      include  'sdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   setvar,palloc
      integer   ilast

      save

c     For T-splines return

      if(tsplfl) return

c     Determine necessary storage for mesh lines and allocate storage

      setvar = palloc(206,'NORMV',  numnp*3,2)

      call setclp(hr(np(43)),ndm,numnp)

      call pnorml(mr(np(32)),mr(np(33)),hr(np(43)),hr(np(206)),
     &            mr(np(78)),nie,ndm,nen,nen1,numnp,numel)


c     For NURBS problems return

      if(nurbfl) return

      setvar = palloc(111,'TEMP1',  numnp,1)

      call pextnda(mr(np(33)),mr(np(111)),ilast)

      setvar = palloc(112,'TEMP2',  max(1,ilast),1)

      call pextndb(mr(np(33)),mr(np(111)),mr(np(112)), ilast)

      call pextndc(mr(np(33)),mr(np(32)),mr(np(111)),mr(np(112)),
     &             mr(np(78)))

      setvar = palloc(112,'TEMP2',   0,1)
      setvar = palloc(111,'TEMP1',   0,1)

      end
