c$Id:$
      subroutine pfullscr(k1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set fullscreen plot

c     Inputs:
c        kl   - Flag to set fullscreen (not used for X11)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'wdata.h'

      logical    wflag
      integer    k1

      save

      if(k1.eq.0) then
        wflag      = .true.
      elseif(wflag) then
      else
        write(*,*) ' Must use UPLOT with zero parameter first'
      endif

      end
