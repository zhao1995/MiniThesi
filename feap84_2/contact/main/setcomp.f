c$Id:$
      subroutine setcomp (npair,cs0,cm0,cp0,hic,csw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set initial penetration flag 'ifpck'             03/05/2007
c       2. Activate contact surface based on time function  14/03/2011
c       3. Add 'csw' to argument of routine                 06/11/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: SET COMmons for a contact Pair

c      Purpose: Update values related to the contactpair in commons

c      Inputs :
c         npair   - # of contact pair
c         cs0(*)  - Contact surfaces control data
c         cm0(*)  - Contact material control data
c         cp0(*)  - Contactpair control data
c         hic(*)  - HIstory Correspondence table

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'prld1.h'
      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_mate.h'
      include  'c_pair.h'
      include  'c_tole.h'

      logical   timfl
      integer   npair,csw, hic((c_lp1+c_lp3),*), kv
      real*8    cs0(nr0,n0c1:nc01,*),cm0(nr0,n0c2:nc02,*)
      real*8    cp0(nr0,n0c3:nc03,*)

      save

      call cdebug0 ('    setcomp',-1)

      timfl = csw.ne.1 .and. csw.ne.14 .and. csw.ne.313

c     Load control flags on common

      ifon   = nint(cp0(3,1,npair))
      iffric = nint(cp0(4,1,npair))
      iftimf = nint(cp0(6,1,npair))
      if(iftimf.gt.0 .and. timfl) then
        if(prldv(iftimf).le.0.0d0) then
          ifon = 0         ! Turns off contact for this pair
        endif
      endif
      ifsolm = nint(cp0(2,2,npair))
      ifdeta = nint(cp0(2,3,npair))
      ifaugm = nint(cp0(2,5,npair))
      ifadhe = nint(cp0(2,7,npair))
      ifintm = nint(cp0(2,10,npair))

c     Load tolerances

      tlipen = cp0(3,6,npair)
      tlopen = cp0(4,6,npair)
      tlouts = cp0(5,6,npair)

c     Set penetration check option

      ifpck  = cp0(3,8,npair).ne.0.0d0

c     Load offsets, set length of ch1,ch2,ch3, # of set, solution driver

      rnpair = nint(cp0(1,-1,npair))
      ofh1   = nint(cp0(2,-1,npair))
      ofh3   = nint(cp0(3,-1,npair))
      lh1    = nint(cp0(4,-1,npair))
      lh3    = nint(cp0(5,-1,npair))
      nset   = nint(cp0(6,-1,npair))
      ndrv   = nint(abs(cp0(1, 0,npair)))

c     Load # of element active till now

      nacte  = nint(cp0(11,-1,npair))

c     Set contact search dimension

      cndm   = nint(cp0(13,-1,npair))

c     Load surfaces information

      nsurf1 = nint(abs(cp0(7,-1,npair)))

      ofs1   = nint(abs(cs0(2,-1,nsurf1)))
      neps1  = nint(abs(cs0(3,-1,nsurf1)))
      dnope1 = nint(abs(cs0(4,-1,nsurf1)))
      ifsty1 = nint(abs(cs0(1, 0,nsurf1)))
      nope1  = nint(abs(cs0(2, 0,nsurf1)))

      nsurf2 = nint(abs(cp0(8,-1,npair)))

      ofs2   = nint(abs(cs0(2,-1,nsurf2)))
      neps2  = nint(abs(cs0(3,-1,nsurf2)))
      dnope2 = nint(abs(cs0(4,-1,nsurf2)))
      ifsty2 = nint(abs(cs0(1, 0,nsurf2)))
      nope2  = nint(abs(cs0(2, 0,nsurf2)))

c     Load material information

      nmat1  = nint(cp0(9,-1,npair))
      if (nmat1.ne.0) then
        ifmty1 = nint(cm0(1,0,nmat1))
        ofm1   = nint(cm0(2,-1,nmat1))
      else
        ifmty1 = 0
        ofm1   = 0
      endif

c     Do the same only if mat # 2 is used

      nmat2  = nint(cp0(10,-1,npair))
      if (nmat2.ne.0) then
        ifmty2 = nint(cm0(1,0,nmat2))
        ofm2   = nint(cm0(2,-1,nmat2))
      else
        ifmty2 = 0
        ofm2   = 0
      endif

c     Load correspondence table for history variables

      do kv = 1,c_lp1
        p1(kv) = hic(kv,npair)
      end do
      do kv = 1,c_lp3
         p3(kv) = hic(c_lp1+kv,npair)
      end do

      end
