c$Id:$
      subroutine pedgex(nspin)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    03/01/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Control routine of data inputs based on edge coordinate

c      Inputs:
c         none      - Data retrieved through common blocks

c      Outputs:
c         none      - Data stored in pointers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'edgdat.h'
      include  'pointer.h'
      include  'print.h'
      include  'sdata.h'
      include  'comblk.h'

      logical   oprt,oprth
      integer   nspin

      save

c     Save print flags

      oprt  = prt
      oprth = prth

c     Set edge spin conditions

      if(espfl) then                  ! Set for 'espi'
        call pespin(hr(np(43)),mr(np(266)),mr(np(190)),
     &              ndm,numnp,nspin,prt,prth)
      endif

      prt  = oprt
      prth = oprth

      end
