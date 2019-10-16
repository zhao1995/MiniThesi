c$Id:$
      subroutine crauto(ns,cs0,ics,lics0)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               R.L. Taylor              October 16, 2001            1.0

c      Acronym: Contact Read AUTO

c      Purpose: Determine number & surface data for boundary slideline

c      Inputs:
c         ns      - First surface number to generate

c      Outputs:
c         ns      - # contact surfaces generated
c         ics(*)  - Connection list for generated surfaces
c         lics0   - Length of array to generate
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_dict.h'
      include  'chdata.h'
      include  'sdata.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   ics(*)
      real*8    cs0(*)

      integer   ns,lics0

      call cdebug0 ('  crauto',-1)

      if(ndm.eq.1) then
        write(*,*) 'NOT IMPLEMENTED'
      elseif(ndm.eq.2) then
        call aslid2c(mr(np(221)),cs0,ics, ns,lics0)
      elseif(ndm.ge.3) then
        call aslid3c(mr(np(221)),cs0,ics, ns,lics0)
      endif

      end
