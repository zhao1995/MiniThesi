c$Id:$
      subroutine crel03 (nsopt,scdat,nope,dnope,ics,emax)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct check on np(161) to .eq. instead of .le. 15/01/2008
c          Write message to 'ilg' and 'iow'
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor           March 23, 1997          1.00

c      Acronym: Contact Read ELEMents

c      Purpose: Input of contact surface element nodes
c                with automatic blending function generation

c      Inputs:
c         nsopt   - # of the Sub-command OPTion
c         scdat   - Sub-Command DATa read
c         nope    - # of NOde Per Element
c         dnope   - Dimension of NOde Per Element

c      Outputs:
c         ics(*)  - Contact element nodal connection array
c         emax    - Element max number found
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'print.h'
      include  'cdata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      integer   nsopt,nope,dnope,ics(*),emax, emaxnew, lip
      real*8    scdat(14)

      save

      call cdebug0 ('      crel03',-1)

c     Set last element number

      emaxnew   = nint(scdat(1))
      if(emaxnew.ne.0) then
        emax = emaxnew
      endif

c     Set pointers for temporary storage. N.B. cannot use palloc here

      lip = 1   + dnope*numel

      if(np(161).eq.0) then
        write(ilg,*) ' NO BLENDING SIDE/NODES: Stop'
        write(iow,*) ' NO BLENDING SIDE/NODES: Stop'
        call plstop()
      endif

c     Input blend segments for 2-d problems

      if(ndm.eq.2) then

c       Set remaining temporary location for scratch arrays

        call crblen(nope,dnope,ics(1),emax, nsopt, scdat(1),
     &              hr(np(161)),mr(np(33)),ics(lip),hr(np(43)))

c     Input blend segments for 3-d problems

      elseif(ndm.eq.3) then

c       Set remaining temporary location for scratch arrays

        if(nope.eq.2) then
          call crblen(nope,dnope,ics(1),emax, nsopt, scdat(1),
     &                hr(np(161)),mr(np(33)),ics(lip),hr(np(43)))
        else
          call crblen3(nope,dnope,ics(1),emax, nsopt, scdat(1),
     &                 hr(np(161)),mr(np(33)),ics(lip),hr(np(43)),
     &                 hr(np(26)),ndm,nen,nen1,numnp,numel,
     &                 prt,prth)
        endif

      endif

      end
