c$Id:$
      subroutine plclos()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Close plot device

c      Inputs:
c         none

c      Outputs:
c         none      - Returns command outputs to text device
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pdata2.h'
      include  'plflag.h'
      include  'print.h'
      include  'x11f.h'

      save

c     Close plot device

      if(.not.fopn) return
      fopn = .false.

c     X11 device

      if(screfl) call gdx11(5,xx,yy)

      end
