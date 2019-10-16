c$Id:$
      subroutine smooth(x,nsmth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Smoothed nodes

c      Inputs:
c         x(ndm,*)  - Nodal coordinates (unsmoothed)
c         nsmth     - Number of smooth cycles

c      Outputs:
c         x(ndm,*)  - Nodal coordinates (smoothed)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none


      include  'cdata.h'
      include  'idptr.h'
      include  'sdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   setvar,palloc
      integer   nsmth, kp, ndfold
      real*8    x(ndm,*)

      save

      setvar = palloc(111,'TEMP1', numnp*ndf, 1)

c     Compute number elements connected to each node

      ndfold = ndf
      ndf    = 1
      call pzeroi(mr(np(111)),numnp)
      call elcnt(numel,nen,nen1,mr(id31),mr(np(33)),mr(np(111)),0)
      call sumcnt(mr(np(111)),numnp,kp)
      ndf    = ndfold

c     Store elements in list

      setvar = palloc(112,'TEMP2', kp, 1)
      call pelcon(numel,nen,nen1,mr(np(33)),mr(id31),
     &            mr(np(111)),mr(np(112)),kp,1)

c     Determine nodes connected to each node

      setvar = palloc(113,'TEMP3', numnp, 1)
      call pnodcn(mr(np(33)),mr(np(78)),mr(np(111)),mr(np(112)),
     &            mr(np(113)),ndm,nen,nen1,numnp,kp)
      setvar = palloc(113,'TEMP3', kp   , 1)

      call xsmth(mr(np(111)),mr(np(113)),x,ndm,numnp,nsmth)

      setvar = palloc(113,'TEMP3', 0, 1)
      setvar = palloc(112,'TEMP2', 0, 1)
      setvar = palloc(111,'TEMP1', 0, 1)
      end
