c$Id:$
      subroutine formgh(ixg,elist)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Forms finite element arrays as required

c      Inputs:
c         ixg(*) - List of elements with slaved nodes
c         elist  - Number of elemenets with slaved nodes

c      Outputs:
c         r(*)       - Export resid.  (TEMP4: hr(np(114)))
c         hh(neqg,*) - Export matrix  (TEMP5: hr(np(115)))
c         g(neq,*,2) - G-vector       (TEMP6: hr(np(116)))
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'eqslv.h'
      include  'eqsym.h'
      include  'hdatam.h'
      include  'idptr.h'
      include  'part0.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   n, nl1, elist, ixg(*)

      save

c     Form appropriate finite element arrays

c             PFORM: ul   ,    xl   ,    tl   ,    ld   ,    p    ,
c                     s   ,    ie   ,     d   ,    id   ,    x    ,
c                    ix   ,    rben ,    ftn  ,     t   ,    jp   ,
c                     u   ,    vel  , dum, dum, dum,

      hflgu  = .false.
      h3flgu = .false.

      call pzero(hr(np(115)),neqg*neqg )
      call pzero(hr(np(116)),neqg*neq*2)
      do n = 1,elist
        nl1 = ixg(n)
        call pform(hr(np(41)),hr(np(44)),hr(np(39)),mr(np(34)),
     &             hr(np(35)),hr(np(36)),mr(np(32)),hr(np(25)),
     &             mr(id31),hr(np(43)),mr(np(33)),mr(np(181)),
     &             hr(np(30)),hr(np(38)),mr(np(20+npart)),hr(np(40)),
     &             hr(np(42)),hr(1),hr(1),hr(1),ndd,nie,ndf,ndm,
     &             nen1,nst,.false.,.false.,.false.,.false.,3,nl1,
     &             nl1,1)

c       Assemble H and R

        fp(1) = np(33) + (ixg(n) - 1)*nen1
        call dasblgh(mr(fp(1)),mr(id31),mr(np(34)),mr(np(112)),
     &               mr(np(113)),hr(np( 35)),hr(np( 36)),
     &               hr(np(114)),hr(np(115)),hr(np(116)))

      end do ! n

      call formhh(hr(np(115)),hr(np(116)),neqg,neq)

      end
