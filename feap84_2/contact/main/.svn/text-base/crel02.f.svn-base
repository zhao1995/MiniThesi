c$Id:$
      subroutine crel02 (nsopt,scdat,nope,dnope,ics,emax,polfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor           March 31, 1996
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact Read ELEMents

c      Purpose: Input of contact surface element nodes
c                with automatic block generation

c      Inputs:
c         nsopt   - # of the Sub-command OPTion
c         scdat   - Sub-Command DATa read
c         nope    - # of NOde Per Element
c         dnope   - Dimension of NOde Per Element
c         polfl   - Flag for polar/cartesian generation

c      Outputs:
c         ics(*)  - Contact element nodal connection array
c         emax    - Element max number found
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'print.h'
      include  'cdata.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      logical   polfl
      integer   nsopt,nope,dnope,ics(*),emax, emaxnew, lip,lep
      real*8    scdat(14)

      save

      call cdebug0 ('      crel02',-1)

c     Set last element number

      emaxnew   = nint(scdat(1))
      if(emaxnew.ne.0) then
        emax = emaxnew
      endif

c     Set pointers for temporary storage. N.B. cannot use palloc here

      lip = 1   + dnope*numel
      lep = lip + numnp

c     Input block segments for 2-d problems

      if(ndm.eq.2) then

c       Set remaining temporary location for scratch arrays

        call crblok(nope,dnope,ics(1),emax, nsopt, scdat(1),polfl,
     &              mr(np(33)),ics(lip),ics(lep),
     &              hr(np(26)),hr(np(43)),ndm,nen,nen1,numnp,numel)

c     Input block segments for 3-d problems

      elseif(ndm.eq.3) then

c       Set remaining temporary location for scratch arrays

        call crblok3(nope,dnope,ics(1),emax, nsopt, scdat(1),polfl,
     &               mr(np(33)),ics(lip),hr(np(43)),
     &               hr(np(26)),ndm,nen,nen1,numnp,numel,prt,prth)

      endif

      end
