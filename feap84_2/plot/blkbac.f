c$Id:$
      subroutine blkbac(blk)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Change PostScript file background to black or white

c      Inputs:
c         blk      - Flag, if true background is black

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   blk

      save

c     Put out commands to make a black background

      if(blk) then
        call fppsin(' 0 0 0 h ')
      else
        call fppsin(' 1 1 1 h ')
      endif
      call fppsdu()
      call fppsin('n 0 0 m 612 0 l 612 792 l 0 792 l c f')
      call fppsdu()
      call fppsin('grestore ')
      call fppsdu()
      if(blk) then
        call fppsin(' -25 0 translate 1 1 1 h ')
      else
        call fppsin(' -25 0 translate 0 0 0 h ')
      endif
      call fppsdu()

      end
