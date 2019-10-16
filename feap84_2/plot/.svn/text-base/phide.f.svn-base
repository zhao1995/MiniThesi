c$Id:$
      subroutine phide(ct,nxd,nxn,nne,nface,iln)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'npix' to 'plix'                          16/11/2011
c       2. Remove arg 3 and arg 4 from call to perspz       09/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Do hidden line removal

c      Inputs:
c         ct        - Plot negative faces if positive
c         iln(2)    - Line type data
c         nface     - Number of faces on surfaces

c      Outputs:
c         nxd       - Face connection dimension
c         nxn       - Number of nodes/face
c         nne       - Number of faces
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comblk.h'
      include  'cdata.h'
      include  'pdata3.h'
      include  'pdatay.h'
      include  'plflag.h'
      include  'pointer.h'
      include  'ppers.h'
      include  'sdata.h'

      integer   nxd,nxn,nne,nface, iln(2)
      real*8    ct

      save

c     Plot visible mesh

      call pzeroi(mr(np(66)),numnp)
      plix    = np(54)
      nxd     = 7
      nxn     = 4
      nne     = nface
      nfac(1) = nne
      call plface(mr(plix),mr(np(62)),hr(np(53)),
     &            3,nxd,numnp,nne,iln,ct)
      if(ndm.eq.3 .and. nen.gt.3) then
        call p3edge(mr(plix),hr(npxx),ndm,nface,nxd)
        edgfl = .true.
      endif

c     Set plot sequence for z-sort

      if(kpers.ne.0) then
        call perspz(hr(np(53)),mr(plix),mr(np(62)),
     7              nxd,nxn,3,numnp,nne)
      endif

      end
