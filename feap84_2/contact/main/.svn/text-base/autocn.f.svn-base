c$Id:$
      subroutine autocn(ix,ib,ip,norm,ma,ndm,nen1,ns)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Auto surface descriptions in 2-d

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include   'c_0.h'
      include   'c_contac.h'
      include   'cdata.h'

      include   'pointer.h'
      include   'comblk.h'

      logical    setvar, palloc
      integer    ma,ndm,nen1,ns

      integer    ix(nen1,numel), ib(numnp), ip(numnp)
      real*8     norm(3,numnp)

      save

      call cdebug0 ('  autocn',-1)

c     Output boundary parameters

      if(ifdb) then
        call iprint(  ib,1,numnp,1,'BOUND?')
        call iprint(  ip,1,numnp,1,'IP')
      endif

c     Build Slideline patches for 2-D problems

      if(ndm.eq.2) then

        setvar = palloc( 221,'ASLD2', 2*numnp  , 1)
        setvar = palloc( 220,'ACON2', 2*numnp  , 1)
        setvar = palloc( 222,'ACIQ2', ip(numnp), 1)

        call aslid2a(mr(np(220)),ib,ip,mr(np(222)),ix,
     &               hr(np(206)),ma,nen,nen1,numnp,numel)

        setvar = palloc( 222,'ACIQ2', 0, 1)
        call aslid2b(mr(np(220)), mr(np(221)),numnp, ns)
        setvar = palloc( 220,'ACON2', 0, 1)

c     Build Slideline patches for 3-D problems

      elseif(ndm.eq.3) then

        call aslid3da(ib,ip,ix,norm,ma,nen,nen1,numnp,numel)

        setvar = palloc( 222,'ACIQ2',     numnp , 1)
        setvar = palloc( 223,'ACON3',2*ip(numnp), 1)
        call aslid3db(mr(np(222)),mr(np(223)),ib,ip,ix,norm,
     &                ma,nen,nen1,numnp,numel,ns)
        setvar = palloc( 223,'ACON3', 0        , 1)
        setvar = palloc( 222,'ACIQ2', 0        , 1)

      endif

      end
