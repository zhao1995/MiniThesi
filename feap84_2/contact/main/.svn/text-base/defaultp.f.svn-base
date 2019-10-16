c$Id:$
      subroutine defaultp (npair,cs0,cm0,cp0,cm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: DEFAULT for contact Pair

c      Purpose: set defaults for contact pair control data

c      Inputs :
c         npair   - # of current pair

c      Outputs:
c         cp0(*)  - Contact pair control data
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_dict.h'
      include  'iofile.h'

      integer   npair, rnsurf1,rnsurf2,rnmat1,rnmat2,kc, ii, opp
      integer   mtyp1,mtyp2, ofsm1,ofsm2
      real*8    cs0(nr0,n0c1:nc01,*),cm0(nr0,n0c2:nc02,*)
      real*8    cp0(nr0,n0c3:nc03,*),cm(*)

      save

      call cdebug0 ('  defaultp',0)

c     On/off of features
c     on normal contact only

      if (nint(cp0(1,1,npair)).eq.0) then
        cp0(3,1,npair) = 1
      elseif (cp0(2,1,npair).eq.2) then
        if (nint(cp0(3,1,npair)).eq.0) then
          cp0(3,1,npair) = 1
        endif
      endif

c     Solution method
c     No default -  penalty term for normal stiffness must be defined

      if (nint(cp0(1,2,npair)).eq.0) then
        write(ilg,3000) (cis(opp(3,2,ii)),ii=1,4)
        write(iow,3000) (cis(opp(3,2,ii)),ii=1,4)
        call plstop()
      endif

c     Detection algorithm
c     Default check all possibilities between two surfaces

      if (nint(cp0(1,3,npair)).eq.0) then
        cp0(2,3,npair) = 1
      endif

c     Material properties
c     Default standard "no material", no friction

      if (nint(cp0(1,4,npair)).eq.0) then
        cp0(3,4,npair) = 1
        cp0(4,4,npair) = 0
      endif

c     Augmentations
c     Default all set to OFF -> 0

      if (nint(cp0(1,5,npair)).eq.0) then
        cp0(2,5,npair) = 1
      else
        if (nint(cp0(2,5,npair)).eq.0) then
          cp0(2,5,npair) = 2
        endif
      endif

c     Tolerances

      if (nint(cp0(1,6,npair)).eq.0) then
        cp0(3,6,npair) = 1.d-8
        cp0(4,6,npair) = 1.d-8
        cp0(5,6,npair) = 1.d-5
      endif

c     Interpolation method: Lagrange

      if (nint(cp0(1,10,npair)).eq.0) then
        cp0(2,10,npair) = 1.0d0
      elseif(nint(cp0(2,10,npair)).eq.0) then
        cp0(2,10,npair) = 1.0d0
      endif

c     Find and store internal number for surfaces and materials

      rnsurf1 = nint(cp0(2,0,npair))
      rnsurf2 = nint(cp0(3,0,npair))
      do kc = 1, cck(1)

c       Check for first surface

        if (cs0(1,-1,kc).eq.rnsurf1) then
          cp0(7,-1,npair) = kc
        endif

c       Check for second surface (allow for single surface contacts)

        if(cs0(1,-1,kc).eq.rnsurf2) then
          cp0(8,-1,npair) = kc
        endif
      end do ! kc

      rnmat1 = nint(cp0(3,4,npair))
      rnmat2 = nint(cp0(4,4,npair))
      do kc = 1, cck(2)

c       Check for first material number

        if (cm0(1,-1,kc).eq.rnmat1) then
          cp0(9,-1,npair) = kc
        endif

c       Check for second material number

        if(cm0(1,-1,kc).eq.rnmat2) then
          cp0(10,-1,npair) = kc
        endif
      end do ! kc

c     If material properties have friction set default to on unless
c     swit off

      if(numcm.gt.0 .and. cp0(2,1,npair).ne.1) then
        rnmat1 = nint(cp0( 9,-1,npair))
        rnmat2 = nint(cp0(10,-1,npair))
        mtyp1  = nint(cm0(1, 0,rnmat1))
        ofsm1  = nint(cm0(2,-1,rnmat1))
        mtyp2  = nint(cm0(1, 0,rnmat2))
        ofsm2  = nint(cm0(2,-1,rnmat2))
        if((mtyp1.gt.0 .and. ofsm1.gt.0) .or.
     &     (mtyp2.gt.0 .and. ofsm2.gt.0)) then
          if(cm(ofsm1).ge.0.0d0 .or. cm(ofsm2).ge.0.0d0) then
            cp0(1,1,npair) = 1
            cp0(2,1,npair) = 2
            cp0(3,1,npair) = 1
            cp0(4,1,npair) = 1
          endif
        endif
      else
        cp0(4,1,npair) = 0
      endif

c     Formats

3000  format(' *ERROR* DEFAULTP: Solution method must be input.'/
     &       '         OPTIONS : ',8(a:,', '))
      end
