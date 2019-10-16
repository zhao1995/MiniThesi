c$Id:$
      subroutine pfullscr(k1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Replace MSFLIB by DFLIB                          08/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set graphics screen to full size or reduced size.

c      Inputs:
c         k1      - Set to full screen if zero
c                 - Set to reduced screen if unity

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      use        DFLIB

      implicit   none

      include   'wdata.h'

      integer         idxl,idyl,jfill
      common /vgraph/ idxl,idyl,jfill

      integer*2       nxpix,nypix,nrows,ncols
      common /vgsize/ nxpix,nypix,nrows,ncols

      logical    wflag
      integer    wxyold(2),xpxold(2) ,wxynew(2),xpxnew(2)
      integer    k1, status, vclrwk

      data       wflag / .false. /

c     Clear the entire screen

      call clearscreen($GCLEARSCREEN)

c     Change to full screen

      if(k1.eq.0) then
        if(.not.wflag) then
          wxyold(1) = wxy(1,2,1)
          wxyold(2) = wxy(2,2,1)
          xpxold(1) = xpxl(1,1)
          xpxold(2) = xpxl(2,1)
          wxynew(1) = nypix*1.47
          wxynew(2) = nypix
          xpxnew(1) = 28*idxl/40
          xpxnew(2) = 28*idyl/40
        endif

        wxy(1,2,1) = wxynew(1)
        wxy(2,2,1) = wxynew(2)

        xpxl(1,1)  = xpxnew(1)
        xpxl(2,1)  = xpxnew(2)
        call settextwindow(int2(nrows-4),int2(1),int2(nrows),int2(65))
        status     = displaycursor($gcursoron)
        status     = setfont('t''Arial''h10b')
        status     = vclrwk()
        wflag      = .true.

c     Change to original screen

      elseif(wflag) then
        status     = vclrwk()
        wxy(1,2,1) = wxyold(1)
        wxy(2,2,1) = wxyold(2)
        xpxl(1,1)  = xpxold(1)
        xpxl(2,1)  = xpxold(2)
        status     = vclrwk()
        call settextwindow(int2((6.6*nrows)/10+1),int2(1),
     &                     int2(nrows),int2(80))
        status     = displaycursor($gcursoron)
        status     = setfont('t''Arial''h10b')

c     User option

      else
        write(*,*) ' Must use UPLOT with zero parameter first'
      endif

      end
