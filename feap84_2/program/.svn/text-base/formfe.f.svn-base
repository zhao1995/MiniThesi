c$Id:$
      subroutine formfe(pnu,pnb,pna,pnl,aufl,bfl,alfl,dfl,
     &                  isw,nl1,nl2,nl3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change np(244) -> np(256) in call to pnforc      02/11/2006
c       2. Remove interupt tests (intf and lintr)           17/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Forms finite element arrays as required

c      Inputs:
c         pnu    - Pointer to current nodal solution vectors
c         pnb    - Pointer for FEM Vector
c         pna    - Pointer for FEM diagonal and upper part of array
c         pnl    - Pointer for FEM lower part of array
c         aufl   - If true assemble 'a' array (which includes 'au')
c         alfl   - If true assemble 'al' array
c         bfl    - If true assemble 'b' array
c         dfl    - If true assembel 'b' uncompressed
c         isw    - Solution switch controlling action to be taken
c         nl1    - First element to be processed
c         nl2    - Last  element to be processed
c         nl3    - Increment to 'nl1'

c      Outputs:
c         hr(pnb)- Values for FEM Vector
c         hr(pna)- Values for FEM diagonal and upper array part
c         hr(pnl)- Values for FEM lower array part
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_tanfl.h'
      include  'cdata.h'
      include  'cdat1.h'
      include  'hdatam.h'
      include  'idptr.h'
      include  'ldata.h'
      include  'ndata.h'
      include  'part0.h'
      include  'pointer.h'
      include  'prflag.h'
      include  'rigid2.h'
      include  'sdata.h'
      include  'comblk.h'

      include  'p_formfe.h'
      include  'p_int.h'

      logical   aufl,bfl,alfl,dfl, flg
      integer   isw, nl1, nl2, nl3

      save

c     Set flag for convergence checks

      if(isw.eq.3 .or. (isw.eq.6 .and. .not.dfl)) then
        lvcn     = lv
        floop(1) = .true.
      endif

c     Copy pointers and flags for contact and rigid bodies

      lafl = alfl
      uafl = aufl
      dbfl = bfl
      ddfl = dfl

      rpu  = pnu
      rpb  = pnb
      rpa  = pna
      rpl  = pnl

c     Form appropriate finite element arrays

c                    ul   ,    xl   ,    tl   ,    ld   ,    p    ,
c                     s   ,    ie   ,     d   ,    id   ,    x    ,
c                    ix   ,    rben ,    ftn  ,     t   ,    jp   ,
c                     u   ,    vel  ,

      call pform(hr(np(41)),hr(np(44)),hr(np(39)),mr(np(34)),hr(np(35)),
     &           hr(np(36)),mr(np(32)),hr(np(25)),mr(id31),hr(np(43)),
     &           mr(np(33)),mr(np(181)),hr(np(30)),hr(np(38)),
     &           mr(np(20+npart)),hr(pnu),hr(np(42)),hr(pnb),hr(pna),
     &           hr(pnl),ndd,nie,ndf,ndm,nen1,nst,aufl,bfl,alfl,dfl,
     &           isw,nl1,nl2,nl3)

c     Add concentrated nodal terms: stiffness, damping, mass

      if(nmfl) then
        call pnodal(mr(id31),hr(np(88)),hr(np(86)),hr(np(87)),
     &              mr(np(29)+ndf*numnp) ,hr(pna),hr(pnb),hr(pnu),
     &              hr(np(42)),ndf*numnp,aufl,bfl,dfl,isw)
      endif

c     Add concentrated radial follower forces

      if(nffl) then
        call pnforc(mr(id31),mr(np(29)),hr(np(256)),
     &              mr(np(20+npart)),hr(pna),hr(pnl),hr(pnb),
     &              hr(np(43)),hr(pnu),alfl,aufl,bfl,dfl,isw)
      endif

c     Interface contributions

      if(np(210).ne.0 .and. isw.ne.8) then

c       Shifts to local arrays for second element

        fp(1) = np(39) + nen       ! TL
        fp(2) = np(41) + ndf*nen*7 ! UL
        fp(3) = np(44) + ndm*nen   ! XL

        call ipform(hr(np(41)),hr(fp(2)),hr(np(44)),hr(fp(3)),
     &              hr(np(39)),hr(fp(1)),mr(np(34)),hr(np(35)),
     &              hr(np(36)),mr(np(32)),hr(np(25)),mr(id31),
     &              hr(np(43)),mr(np(33)),mr(np(210)),hr(np(30)),
     &              hr(np(38)),mr(np(20+npart)),hr(pnu),hr(np(42)),
     &              hr(pnb),hr(pna),hr(pnl),aufl,bfl,alfl,dfl,isw)
      endif

c     Form rigid body arrays

      call rigidb (2,isw,flg)

c     Contact contributions

      call contact (isw)

c     User contributions

      call uformfe(pnu,pna,pnl,pnb,
     &             aufl,alfl,bfl,dfl,isw,nl1,nl2,nl3)

c     Peform final step for parallel solutions

      call parform(isw)

c     Reset update flag for history variables

      hflgu  = .false.
      h3flgu = .false.

      end
