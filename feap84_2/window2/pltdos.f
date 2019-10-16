c$Id:$

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Replace MSFLIB by DFLIB                          08/12/2006
c       2. Definei array 'mask' for 'vsltyp'                01/05/2007
c-----[----------------------------------------------------------------]
c      Purpose:  Library of plot outputs for Windows systems

c      Inputs:
c         See individual routines

c      Outputs:
c         See individual routines
c-----[----------------------------------------------------------------]

      integer function vopnwk()

      use DFLIB

      implicit  none

      include  'iofile.h'
      include  'plflag.h'
      include  'pdata2.h'

      integer         idxl,idyl,jfill
      common /vgraph/ idxl,idyl,jfill

      type(windowconfig) :: myscreen
      type(qwinfo)       :: winfo,frmwindow
      integer(2)         :: numfonts,status
      integer            :: uchild

c     Open workstation, home cursor, set up scaling

      if(screfl) then

        uchild = 20
        open(unit=uchild,file='user',title='F E A P     P l o t s')

c       Get window size data for plot outputs

        status  = getwsizeqq(uchild,qwin$sizecurr,winfo)
        status  = getwsizeqq(qwin$framewindow,qwin$sizemax,frmwindow)

c       Position for Graphics screen (in text row/columns shifts)

        winfo.X = frmwindow.X  ! frmwindow.W - winfo.W
        winfo.Y = frmwindow.Y  ! frmwindow.H - winfo.H

        winfo.type=qwin$set

        status  = setwsizeqq(uchild,winfo)
        if(status.ne.0) status = setwsizeqq(uchild,winfo)

c       Get screen capability

        status = getwindowconfig(myscreen)

c       Set sizes to maximum available

        myscreen.numtextcols=-1
        myscreen.numtextrows=-1
        myscreen.fontsize   =-1

c       Sizing graphics screen

        myscreen.numypixels= (myscreen.numypixels)*0.6
        myscreen.numxpixels= (myscreen.numypixels)*1.3

c       Set window configuration

        status = setwindowconfig(myscreen)
        if(status.ne.0) status = setwindowconfig(myscreen)

c       Set sizing for lines drawn by FEAP

        idyl   = nint(22480.0/(myscreen.numypixels))
        idxl   = idyl

        if(myscreen.numcolors .le. 4) then
          jfill = 1
        else
          jfill = 2
        endif

        vopnwk = displaycursor ( $GCURSOROFF   )
        call     clearscreen   ( $GCLEARSCREEN )

c       Set font for Arial outputs vector mode

        numfonts=initializefonts()
        if (numfonts.le.0) print *,"INITIALIZEFONTS error"
        if (grstatus().ne.$GROK) then
          write(*,*) 'INITIALIZEFONTS GRSTATUS error.'
        endif
        vopnwk =  setfont( 't''Arial''h14b' )

      endif

c     Tile Windows for Maximum Viewing

      status = clickmenuqq(QWIN$TILE)
      status = focusqq(0)
      status = setactiveqq(uchild)

      vopnwk = 0

      end

      integer function vtxsiz(isw)

c     Set font size for Helvetica Bold outputs vector mode

      use DFLIB

      implicit  none

      integer   isw

      save

      if(isw.eq.1) then
        vtxsiz =  setfont( 't''Arial''h14b' )
      elseif(isw.eq.2) then
        vtxsiz =  setfont( 't''Arial''h16b' )
      elseif(isw.eq.3) then
        vtxsiz =  setfont( 't''Arial''h20b' )
      endif
      vtxsiz = 0

      end

      integer function vtxwin()

      vtxwin = 0

      end

      integer function vclrwk()

c     Clear the workstation

      use DFLIB

      implicit  none

      save

      call clearscreen( $GCLEARSCREEN )
      vclrwk = 0

      end

      integer function vhomwk()

c     Home cursor - text mode

      use DFLIB

      implicit  none

      type(rccoord) s

      save

      call settextposition( int2(1) , int2(1) , s )
      vhomwk = 0

      end

      integer function vclswk()

c     Function to close plot (if necessary)

      implicit  none

      save

      close(20, status='delete')

      vclswk = 0
      end

      integer function vgtxts(xi,yi,nn,cstr)

c     Place graphics text on screen

      use DFLIB

      implicit  none

      integer         idxl,idyl,jfill
      common /vgraph/ idxl,idyl,jfill

      integer          :: n,nn
      integer(2)       :: ix,iy
      real(8)          :: xi,yi
      character(len=1) :: cstr(nn)
      type(xycoord)    :: xy

      save

c     x,y locations for outgtext

      ix = xi*22000/idxl
      iy = (22200 - yi*22000)/idyl

      call moveto( ix , iy , xy )

c     Output characters one at a time

      do n = 1,nn
         call outgtext(cstr(n))
      end do ! n

      vgtxts = 0

      end

      integer function vsltyp(it)

      use DFLIB

      implicit  none

      integer     it
      integer*2   mask(7)

      save

c     Set line patterns (1-bit draws; 0-skips)

      data        mask /  2#1111111111111111 ,
     &                    2#1111111100000000 ,
     &                    2#1100110011001100 ,
     &                    2#1111000011110000 ,
     &                    2#1000100010001000 ,
     &                    2#1111001100110010 ,
     &                    2#1111100011111000 /

      call setlinestyle( mask(it) )
      vsltyp = 0

      end

      integer function vipal(it)

      use       DFLIB

      implicit  none

      integer         idxl,idyl,jfill
      common /vgraph/ idxl,idyl,jfill

      integer :: it

      integer    ipal(18)

      save

      data  ipal/  #FFFFFF      ,   !   1: BRIGHTWHITE
     &             #0000FF      ,   !   2: RED
     &             #00FF00      ,   !   3: GREEN
     &             #FF0000      ,   !   4: BLUE
     &             #00FFFF      ,   !   5: YELLOW
     &             #FFFF00      ,   !   6: CYAN
     &             #FF00FF      ,   !   7: MAGENTA
     &             #00A5FF      ,   !   8: ORANGE
     &             #507FFF      ,   !   9: CORAL
     &             #2FFFAD      ,   !  10: GREEN YELLOW
     &             #B3DEF5      ,   !  11: WHEAT
     &             #FF6941      ,   !  12: ROYAL BLUE
     &             #E020A0      ,   !  13: PURPLE
     &             #D4FF7F      ,   !  14: AQUAMARINE
     &             #9020D0      ,   !  15: VIOLET RED
     &             #8B3D48      ,   !  16: DARK SLATE BLUE
     &             #BEBEBE      ,   !  17: GRAY
     &             #D3D3D3      /   !  18: LIGHTGRAY

c     Set color pallet

      if(it.gt.0 .and. it.le.18 ) then
        vipal = ipal(it)
        if(jfill.lt.2) vipal = 1
      else
        vipal = #000000 ! Black
      endif

      end

      integer function vstcol(it)

      use DFLIB

      implicit  none

c     Set text color for graphics output

      integer    :: icll
      integer    :: it, vipal

      save

      icll   = vipal(it)
      vstcol = settextcolorrgb( icll )

      end

      integer function vslcol(it)

      use DFLIB

      implicit  none

c     Set line color for graphics output

      integer    :: icll
      integer    :: it, vipal

      save

      icll   = vipal(it)
      vslcol = setcolorrgb( icll )

      end

      integer function vpline(ixy,ipen)

c     Move/draw for lines

      use DFLIB

      implicit  none

      integer         idxl,idyl,jfill
      common /vgraph/ idxl,idyl,jfill

      type(xycoord) :: xy
      integer(2)    :: ix,iy
      integer       :: ipen
      integer       :: ixy(2,*)

      save

c     Set cocordinates

      ix = ixy(1,1)/idxl
      iy = (22000 - ixy(2,1))/idyl

c     Draw line

      if(ipen.eq.2) then
        vpline = lineto( ix , iy )

c     Move without draw

      elseif(ipen.eq.3) then
        call     moveto( ix , iy , xy )
      end if

      end

      integer function vfarea(npt,ixy)

      use DFLIB

c     Panel fill

      implicit  none

      integer         idxl,idyl,jfill
      common /vgraph/ idxl,idyl,jfill

      type(xycoord) :: poly (62)
      integer       :: npt, ixy(2,npt)
      integer(2)    :: n, nn

      save

c     Trace area to fill

      nn = min(31,npt)
      do n = 1,nn
        poly(n).xcoord = ixy(1,n)/idxl
        poly(n).ycoord = (22000 - ixy(2,n))/idyl
      end do ! n

c     Perform fill

      vfarea = polygon( $GFILLINTERIOR, poly , nn )
      vfarea = 0

      end
