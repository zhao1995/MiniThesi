c----------------------------------------------------------------------
c.... Additional programs for WIN-INTEL
c     and dummy-subroutines for UNIX,PHIGS,GKS,WIN_SALFORD
c
c-----------------------------------------------------------------------
c.... INTEL Graphik - Routinen


c---- SALFORD ----------------------------------------------------------
c
c     SUPERLUINT
c
      subroutine datri3(a,ia,neq,nneg)
      end
c
      subroutine dasol3 (a,b,ia,neq,energy)
      end
c


c---- INTEL FINAL---------------------------------------------------------
c
c-----------------------------------------------------------------------
      subroutine clear_screen
c-----------------------------------------------------------------------
c
c.... Purpose: clear screen
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      call clearscreen ($GCLEARSCREEN)
#endif
      return
      end
c
c-----------------------------------------------------------------------
      subroutine clear_screen_area(ix1,iy1,ix2,iy2,icol)
c-----------------------------------------------------------------------
c
c.... Purpose: clear screen area
c
c     Input:
c     ix1,iy1   - LL corner rectangle
c     ix2,iy2   - UR corner rectangle
c     icol      - color
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      integer*2 ix1,iy1,ix2,iy2,icol
cww intel clear all:  open, at the moment:
      call clearscreen ($GCLEARSCREEN)
#endif
      return
      end
c
c-----------------------------------------------------------------------
      subroutine draw_line(ix1,iy1,ix2,iy2,icol)
c-----------------------------------------------------------------------
c
c.... Purpose: draw line from point 1 to point 2
c
c     Input:
c       ix1,iy1   - point 1
c       ix2,iy2   - point 2
c       icol      - color
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      integer*2 ix1,iy1,ix2,iy2,icol
      integer*4 result
      type (xycoord) xy
cww   icol=???  Farbe setzen!!
      call moveto(ix1,iy1,xy)
      result = lineto(ix2,iy2)
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine draw_text(str,ixa,iya,icc)
c-----------------------------------------------------------------------
c
c.... Purpose: draw text
c
c     Input:
c     str      - string
c     ixa,iya  - position
c     icc      - color
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c
#ifndef _BATCH_
      USE ifqwin
      integer*2 ixa,iya,icc
      character*(*) str
      integer*4 result
      type (xycoord) xy
cww   icol=???  Frabe setzen!!
      call moveto(ixa,iya,xy)
      Call outgtext (str)
#endif
      end
c
c-----------------------------------------------------------------------
c
      logical(4) function initialsettings( )
c-----------------------------------------------------------------------
c
c.... Purpose: initial settings for FEAP-Intel
c              set frame to maximum size
c              set about box
c
c     Input:
c
c     Output:
c       frame, aboutbox
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c
#ifndef _BATCH_
      USE ifqwin
      logical(4) resultl
      logical statusl
      integer(4) result
      type(qwinfo) winfo
      external mesh,macro,plot,elem,nege1,nege2,gmesh1,gmesh2,cadobe

c.... Set the About box message
      result = aboutboxqq (
     + ' F E A P\r  Finite Element Analysis Program\r
     +  (c) Institut für Baustatik\r
     +  Karlsruhe Institute of Technology (KIT)\r
     +      Version 01.15'C)


c.... set frame to maximum size
      winfo.type = qwin$max
      statusl    = setwsizeqq(qwin$framewindow,winfo)

c.... set own menue here
c.... extend menue
      resultl=insertmenuqq(7,0,$menuenabled, 'Help-Manual'c, nul)
      resultl=appendmenuqq(7,  $menuenabled, 'Mesh 'c, mesh)
      resultl=appendmenuqq(7,  $menuenabled, 'Macro'c, macro)
      resultl=appendmenuqq(7,  $menuenabled, 'Plot 'c, plot)
      resultl=appendmenuqq(7,  $menuenabled, 'Elem'c, elem)
      resultl=appendmenuqq(7,  $menuseparator, ''c, nul)
      resultl=appendmenuqq(7,  $menuenabled, 'Nege-Input'c, nege1)
      resultl=appendmenuqq(7,  $menuenabled, 'Nege-Examples'c, nege2)
      resultl=appendmenuqq(7,  $menuseparator, ''c, nul)
      resultl=appendmenuqq(7,  $menuenabled, 'Gmesh-Input'c, gmesh1)
      resultl=appendmenuqq(7,  $menuenabled, 'Gmesh-Examples'c, gmesh2)
      resultl=appendmenuqq(7,  $menuseparator, ''c, nul)
      resultl=appendmenuqq(7,  $menuenabled, 'PDF-Manual'c, cadobe)


      initialsettings= .true.
#endif
      end function initialsettings
c
      subroutine  mesh()
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      character name*20
      name = 'mesh'
      call start_manual(name)
#endif
      return
      end
c
      subroutine  macro()
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      character name*20
      name = 'macro'
      call start_manual(name)
#endif
      return
      end
c
      subroutine  plot()
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      character name*20
      name = 'plot'
      call start_manual(name)
#endif
      return
      end
c
      subroutine  elem()
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      character name*20
      name = 'elmt'
      call start_manual(name)
#endif
      return
      end
c
      subroutine  nege1()
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      character name*20
      name = 'nege'
      call start_manual(name)
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine  nege2()
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      character name*20
      name = 'nege_bsp'
      call start_manual(name)
#endif
      return
      end
c
      subroutine  gmesh1()
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      character name*20
      name = 'gmesh'
      call start_manual(name)
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine  gmesh2()
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      character name*20
      name = 'gmesh_bsp'
      call start_manual(name)
#endif
      return
      end
c
      subroutine  cadobe()
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      character name*20
      name = 'titl'
      call start_manual(name)
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pwopn
c-----------------------------------------------------------------------
c
c.... Purpose: open Data window
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      USE wincfg

      implicit none
      type(qwinfo) winfo
      type(windowconfig) feapdialog
      integer*4 status,if1,ipos
      integer*2 nxspix,nyspix,nrows,ncols
      integer*2 nx1pix,ny1pix
      integer(2) numfonts
      character*80 titlei
      logical statusl

c     Open Child Input Window
      open (unit=0,file='con')

c.... Get Screen Size
      status = getwindowconfig(feapdialog)
      nxspix = feapdialog.numxpixels
      nyspix = feapdialog.numypixels
      ncols  = feapdialog.numtextcols
      nrows  = feapdialog.numtextrows

c.... Set Input Data Window
c.... Position and Size
      winfo.type = QWIN$SET
      winfo.x    = nint(fpx1 * dble(ncols))
      winfo.y    = nint(fpy1 * dble(nrows))
      winfo.w    = nint(fwx1 * dble(ncols))-3 ! -3 to have scrollbars
      winfo.h    = nint(fwy1 * dble(nrows))-7
      status   = setwsizeqq(0,winfo)

      feapdialog.fontsize     = #0008000E  ! klein ww
c     feapdialog.fontsize     =  524304   ! groesser mk  02.06.2010

      if1=ipos(title1,40)
      titlei(1:if1) = title1(1:if1)
      titlei(if1+1:80)=' '
      feapdialog.title= titlei
      feapdialog.title= titlei
c     feapdialog.title        = "F E A P  D i a l o g  W i n d o w"C
      feapdialog.bitsperpixel = -1
      feapdialog.mode         = QWIN$SCROLLDOWN

      status = setwindowconfig(feapdialog)
      if (status.eq.0) status = setwindowconfig(feapdialog)

      statusl = displaycursor($gcursoron)

c.... background and text color - simple version
      status = SETBKCOLOR(15)  ! 15=bright white,  original values
      status = SETTEXTCOLOR(0) ! black,            original values
      CALL CLEARSCREEN ($GCLEARSCREEN)

c.... Set font for Arial Bold outputs vector mode
      numfonts=initializefonts()
      if (numfonts.le.0) print *,"INITIALIZEFONTS error"
c     if (grstatus().ne.$GROK) print *,'INITIALIZEFONTS GRSTATUS error.'
      status = setfont( 't''Arial''h16b' )

c.... set help window
cww   call pwopnh

c.... set focus on dialog window
      status =  FOCUSQQ(0)
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pwopnh
c-----------------------------------------------------------------------
c
c.... Purpose: open HELP window
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      USE wincfg

      implicit none
      type(qwinfo) winfo
      type(windowconfig) feaphelp
      integer*4 status
      integer*2 nxspix,nyspix,nrows,ncols

c     Open Child Help Window
      open(2,file="USER", TITLE= 'F E A P  H E L P  W I N D O W')

c.... Get Screen Size
      status = getwindowconfig(feaphelp)
      nxspix = feaphelp.numxpixels
      nyspix = feaphelp.numypixels
      ncols  = feaphelp.numtextcols
      nrows  = feaphelp.numtextrows

c.... Set Help Window
c.... Position and Size
      winfo.type = QWIN$SET
      winfo.x    = nint(fpx3 * dble(ncols))
      winfo.y    = nint(fpy3 * dble(nrows))
      winfo.w    = nint(fwx3 * dble(ncols))-3 ! -3 to have scrollbars
      winfo.h    = nint(fwy3 * dble(nrows))-7
      status   = setwsizeqq(2,winfo)

c.... background and text color - simple version
      status = SETBKCOLOR(15)  ! 15=bright white,  original values
      status = SETTEXTCOLOR(0) ! black,            original values
      CALL CLEARSCREEN ($GCLEARSCREEN)
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plotstwin
c-----------------------------------------------------------------------
c
c.... Purpose: open control window
c
c     Comment: in SR PWOPN
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plotgrwin(isw)
c-----------------------------------------------------------------------
c
c.... Purpose: open graphic window
c
c     Inputs:
c       1  - open  window
c       2  - close window
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      USE iwinio
      USE wincfg
      implicit none

      type(qwinfo) winfo
      type(windowconfig) feapgraph
      integer isw,if1,ipos
      integer status
      integer*2 nrows,ncols
      integer(2) numfonts
      character*80 title2i

      if(isw.eq.2) return

c     Open Child Graphic Window
      open (unit=1,file='USER')

c.... Get Screen Size
      status = getwindowconfig(feapgraph)
      iwxs = feapgraph.numxpixels
      iwys = feapgraph.numypixels
      ncols  = feapgraph.numtextcols
      nrows  = feapgraph.numtextrows

c.... real pixels of graphic window
      iwxgs = iwxs*fwx2
      iwygs = iwys*fwy2

c.... Set Input Graphic Window
c.... Position and size
      winfo.type = QWIN$SET
      winfo.x    = nint(fpx2 * dble(ncols))
      winfo.y    = nint(fpy2 * dble(nrows))
      winfo.w    = nint(fwx2 * dble(ncols))
      winfo.h    = nint(fwy2 * dble(nrows))
      status   = setwsizeqq(1,winfo)

      if1=ipos(title2,40)
      title2i(1:if1) = title2(1:if1)
      title2i(if1+1:80)=' '
      feapgraph.title = title2i
      feapgraph.title = "F E A P  G r a p h i c   W i n d o w"C

      status = setwindowconfig(feapgraph)
      if (status.eq.0) status = setwindowconfig(feapgraph)

c.... background color
      status = SETBKCOLOR(0)  ! original black for LOGO FEAP
      CALL CLEARSCREEN ($GCLEARSCREEN)

c.... Set font for Arial Bold outputs vector mode
      numfonts=initializefonts()
      if (numfonts.le.0) print *,"INITIALIZEFONTS error"
c      if (grstatus().ne.$GROK) print *,'INITIALIZEFONTS GRSTATUS error.'
      status = setfont( 't''Arial''h16b' )
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
cww        subroutine wait_feap(ityp)
c-----------------------------------------------------------------------
c
c.....Purpose:  wait and possible controlled abort, if SR start_process
c               does not answer
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
cww      character*1 y
cww
cww      write(*,*) 'F E A P  generates a mesh.'
cww      write(*,*) 'Wait until the generation is finished by an ouput ',
cww     +           'of generated nodes and elements!'
cww      write(*,*) 'Continue calculation? (y or n)'
cww       read(*,1000) y
cww      if(y.eq.'y'.or.y.eq.' ') then
cww        ityp = i
cww        return
cww      else
cww        stop 'SR wait_feap'
cww      end if
cww
cww1000  format(a1)
cww      end
c
c----------------------------------------------------------------------
c
       integer function license()
c-----------------------------------------------------------------------
c....  Purpose:  Expiration date of FEAP license reached
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c
      write(*,*) 'Expiration date of FEAP-License reached.'
      write(*,*) 'Please ask for a new FEAP_Registration_Key at'
      write(*,*) 'Institut für Baustatik'
      write(*,*) 'Karlsruhe Institute of Technology (KIT)'
      write(*,*) 'Email: info@ibs.kit.edu'
      license = 1
      end
c
c----------------------------------------------------------------------
c
      subroutine center_mouse(deltax,deltay,idev)
c----------------------------------------------------------------------
c
c.....Purpose: Define center of window for Graphical output
c
c     Inputs:
c       deltax  -  0...1
c       deltay  -  0...1
c
c...  Comment: approximately correct
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      USE ifwin
      USE iwinio
      real*8 deltax,deltay
      integer*4 mouseevent,keystate, result
      integer*4 cursor,newcursor
      integer*4 ixm,iym
      if(idev.ne.3) return

c.... set start position of cursor
      ixm = 0
      iym = 0

      mouseevent = MOUSE$LBUTTONDOWN

c...  set mouse cursor to cross hair
      newcursor = LoadCursor(0_HANDLE, IDC_CROSS)
      cursor    = SetMouseCursor(newcursor)

      write(*,*) 'LEFT Button: Define new center of Figure'

      result = WAITONMOUSEEVENT (mouseevent,keystate,ixm,iym)
      deltax =      ixm/(iwxgs*0.75)
      deltay = 1. - iym/(iwygs*1.00)

c...  set mouse cursor to standard arrow
      newcursor = LoadCursor(0_HANDLE, IDC_ARROW)
      cursor    = SetMouseCursor(newcursor)
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sleepw(t)
c-----------------------------------------------------------------------
c
c.... Purpose: Suspends execution of a process for a specified interval
c
c     Inputs
c       t    time in seconds to wait
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE ifport
      real*4 t
      integer*4 hold_time

      hold_time = t
      call sleep(hold_time)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine mess_win(text,ityp)
c----------------------------------------------------------------------
c
c     Purpose: draw message on window
c
c     Inputs:
c     text           message to plot
c     ityp           window  action
c     ityp = -3      text    wait until click                          |
c     ityp = -2      text    wait some times                           |
c     ityp = -1 or 0 text    wait until click error                    |
c     ityp =  1      graphic wait until click error                    |
c     ityp =  2      graphic wait some times                           |
c     ityp =  3      graphic wait until click                          |
c     ityp =  4      open  graphic                                     |
c     ityp =  5      close graphiC                                     |
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      USE ifport
      character*(*) text
      character*229 caption
      integer*4 ityp,result,mtype,hold_time
      hold_time = 1

      if(ityp.eq.-1.or.ityp.eq.0.or.ityp.eq.1) then
        caption = 'F E A P   Problem Message'
        mtype = MB$ICONEXCLAMATION
        result = MESSAGEBOXQQ (text,caption,mtype)
      else if(ityp.eq.-2.or.ityp.eq.2) then
        caption = 'F E A P   Information'
        mtype = MB$ICONEXCLAMATION
cww     result = MESSAGEBOXQQ (text,caption,mtype) ! geht nur mit click zu
        write(*,*) text ! provisorial
      else if(ityp.eq.-3.or.ityp.eq.3) then
        caption = 'F E A P   Message'
        mtype = MB$ICONEXCLAMATION
        result = MESSAGEBOXQQ (text,caption,mtype)
      else if(ityp.eq.4) then
        caption = 'F E A P   Input processing....'
      else if(ityp.eq.5) then

      endif
#endif
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine filnamw
c----------------------------------------------------------------------
c
c.....Purpose: interactive file handling control for FEAP
c              for Windows/INTEL
c
c     Inputs: -
c
c     Outputs: File feapname with 4 entries for inp/outp/rest/plot
c
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      USE comfil
      USE feapprog
      USE iodata
      USE psize
      USE vdata
      integer         ip,ipos
      character*229   fcis
      character*229   feapname

c.... read filenames via external program

      ip=ipos(fpath,229)
      fcis = ' '
      fcis(1:ip) = fpath(1:ip)
      fcis(ip+1:ip+13) = 'filnamw_i.exe'

      call start_process(fcis)

c.... set filenames in FEAP

      feapname= './feapname'
      open(ios,file=feapname,status='unknown')
      read(ios,1000) finp,fout,fres,fsav,fplt
      close(ios)

c.... temporary control of files
      write(*,2001) finp,fout,fres,fsav,fplt

c.... write head
      write(*,2000) versn,maxm

1000  format(5(a229,/))
2000  format(///
     1 '    F I N I T E   E L E M E N T   A N A L Y S I S',
     2        '   P R O G R A M'/
     3        12x,'(C) R.L. Taylor, University of California 2015'
     4              //23x,'VERSION: ',a16/32x,a16/32x,a16/
     5                23x,'STORAGE: ',i12)
2001  format(/1x,'Files are set as:  Filename'/
     1        1x,'Input   (read ) : ',1x,a80/
     2        1x,'Output  (write) : ',1x,a80/
     3        1x,'Restart (read ) : ',1x,a80/
     4        1x,'Restart (write) : ',1x,a80/
     5        1x,'Plots   (write) : ',1x,a80)

      return
      end
c
c----------------------------------------------------------------------
c
      subroutine fexit()
c-----------------------------------------------------------------------
c
c.... Purpose: skip last exit window
c
c     W. Wagner BS KIT 05/10
c-----------------------------------------------------------------------
c
#ifndef _BATCH_
      USE ifqwin                                  !mk  26.05.2010
      result1 = SETEXITQQ (QWIN$EXITNOPERSIST)    !mk 26.05.2010
#endif
      end
c
c-----------------------------------------------------------------------
c
cww   subroutine clwopen(title,ix,iy,ixw,iyw,iw1,iw2)
      subroutine clwopen(title,iw1,iw2)
c-----------------------------------------------------------------------
c
c.... Purpose: open help text window
c
c     W. Wagner BS KIT 07/10
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      integer*4 status
c.... set focus on help window
      status =  FOCUSQQ(2)
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine clwclose(iw1,iw2,istop)
c-----------------------------------------------------------------------
c
c.... Purpose: close input text window
c
c     W. Wagner BS KIT 07/10
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE ifqwin
      integer*4 status
c.... set focus on dialog window back
      status =  FOCUSQQ(0)
#endif
      return
      end
c
c-----------------------------------------------------------------------

c###TO DO ##################################

c-----------------------------------------------------------------------
c
      subroutine copyclip(k1)
c-----------------------------------------------------------------------
c
c.... Purpose: copy graphic window to clipboard
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
cww      include <windows.ins>
#ifndef _BATCH_
      USE iwinio
      integer*4 icopy,x1,y1,x2,y2,k1
      integer*2 ix1,iy1,ix2,iy2
      integer*2 ihg,ivg,bstat
c.... initial values
      x1 = -5
      y1 = 5
      x2 = iwxgs-15
      y2 = iwygs+5
      ihg=0
      ivg=0
      ix1=x1
      iy1=y1
      ix2=x2
      iy2=y2
c
      if(k1.eq.1) then
c....   copy whole window
      elseif(k1.eq.2) then
c....   copy part of window defined by mouse
cww     call clwopen('COPY   Parameters',1,iwys-120,680,150,1,2)
        call clwopen('COPY   Parameters',1,2)
        write(*,2002)
c....   set rubberband
        call set_graphics_selection(1)
c....   find point 1
140     call get_mouse_position(ihg,ivg,bstat)
        if(bstat.eq.0) goto 140 ! no button depressed
        x1 = ihg
        y1 = ivg
        ix1 = x1
        iy1 = y1
c.....  find point 2
150     call get_mouse_position(ihg,ivg,bstat)
        if(bstat.eq.1) goto 150 ! left button depressed
        call get_mouse_position(ihg,ivg,bstat)
        x2 = ihg
        y2 = ivg
        ix2 = x2
        iy2 = y2
c....   set rubberband back
        call set_graphics_selection(0)
        call clwclose(1,2,0)
      endif
c.... plot actual rectangle (like rubberband, only control)
cww      call draw_line@(ix1,iy1,ix2,iy1,6)
cww      call draw_line@(ix2,iy1,ix2,iy2,6)
cww      call draw_line@(ix2,iy2,ix1,iy2,6)
cww      call draw_line@(ix1,iy2,ix1,iy1,6)
cww      icopy = graphics_to_clipboard@(x1,y1,x2,y2)
      write(*,2001)
      return
2001  format(' FEAP-Window is copied to Clipboard')
2002  format(' Mark 1. copy point with left mouse button,',
     +       ' hold button depressed and move ',//,
     +       ' mouse cursor to 2. copy point.',
     +       ' Loose mouse button for chosen rectangle.')
#endif
      end
c
c-----------------------------------------------------------------------
c
      subroutine efeap
c-----------------------------------------------------------------------
c
c     Purpose:  internal procedure editor under WIN/INTEL
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE feapprog

      call start_process(fproce)

      return
      end
c

      integer*4 function open_nf()
      open_nf=1
      end
c
      integer*4 function open_f()
      open_f=1
      end
c
       integer function aboutp()
c----------------------------------------------------------------------
c....  function for  ABOUTP menu                                      |
c----------------------------------------------------------------------
      aboutp = 1
      end
c----------------------------------------------------------------------
c
      integer*4 function helpp()
c----------------------------------------------------------------------
c.....FEAP Procedure Editor Short description
c----------------------------------------------------------------------
      helpp = 1
      return
      end
c
      subroutine opencu(ios,path)
c.....open file
      integer*4 ios
      character*229 path
      open(ios,file=path,status='unknown',form='unformatted')
      return
      end
c
      subroutine opencf(ios,path)
c.....open file
      integer*4 ios
      character*229 path
      open(ios,file=path,status='unknown',form='formatted')
      return
      end
c----------------------------------------------------------------------
c
      subroutine perform(i1,i2,i3,i4)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine perform1(isw,name)
c----------------------------------------------------------------------
c.... show actual macro in data  input
c     WW BS UKA 2/04
c----------------------------------------------------------------------
      return
      end
c
      subroutine perform2(title,isw,nit,nmax)
c----------------------------------------------------------------------
c.... show actual state of  data input for coor and elem
c.....isw=1 open isw=2 show isw=3 close
c.....title: character*4
c     WW BS UKA 03/04
c----------------------------------------------------------------------
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine prin_ps()
c-----------------------------------------------------------------------
c.... print existing postscript or hpgl files
c-----------------------------------------------------------------------
      character*229   newf
cww      integer winio@,i
      external istatps
      newf=' '
c.... look for files
cww      i=winio@('%bg[grey]&')
cww      i=winio@('%ca[Existing PS or HPGL - files in Problem]&')
cww      i=winio@('%ob[depressed]%nl%`ft@@&',' PS-files','*.eps')
cww      i=winio@('%ft@@&','All   files','*.*')
cww      i=winio@('%`^12tt[&PS]&','file_openr',newf,istatps)
cww      i=winio@('%`ft@@&','PGL-files','*.pgl')
cww      i=winio@('%ft@@&','All   files','*.*')
cww      i=winio@('%^12tt[ &HPGL ]&','file_openr',newf,istatps)
cww      i=winio@('%12tt[ &EXIT ]%cb')
      return
      end
c
      integer function  istatps()
      istatps=1
      end
c
      subroutine file_rest(newf)
c-----------------------------------------------------------------------
c.... print existing rest files
c-----------------------------------------------------------------------
      character*229   newf
cww      integer winio@,i
      external irest
      newf=' '
c.... look for files
cww      i=winio@('%bg[grey]&')
cww      i=winio@('%ca[Problem]&')
cww      i=winio@('%ob[depressed]%nl%`ft@@&',' Restart-Files','r*')
cww      i=winio@('%ft@@&','All   files','*.*')
cww      i=winio@('Restart File does not exist!%nl%nl&')
cww      i=winio@('%`^8tt[&Files]%cb','file_openr',newf,irest)
      return
      end
c
      integer function  irest()
      irest=-1
      end
c
c----------------------------------------------------------------------
c
      subroutine readstr(text1,text2)
c-----------------------------------------------------------------------
c.....  input character string                                       |
c-----------------------------------------------------------------------
      return
      end
c
c-----------------------------------------------------------------------
c
      integer function  istati()
c-----------------------------------------------------------------------
c
c.... test if input file exists
c-----------------------------------------------------------------------
c
      USE comfil
      USE istat
      logical linp,pcomp
      external istato,istatr,istatw
      if(newf.ne.' ') then
        finp = newf
        newf = ' '
      endif
      wdi = 'new'
      inquire(file=finp,exist=linp)
      if(linp) wdi = 'exists'
c.... set default files for a filname beginning with 'I'
      call dochar2(finp,ipos)
      if(pcomp(finp(ipos:ipos),'I',1)) then
        fout = finp
        fres = finp
        fsav = finp
        fplt = finp
        call dochar1(fout,'O',ipos)
        call dochar1(fres,'R',ipos)
        call dochar1(fsav,'R',ipos)
        call dochar1(fplt,'P',ipos)
        i1 = istato()
        i1 = istatr()
        i1 = istatw()
      endif

      istati=1
      end
c
c-----------------------------------------------------------------------
c
      integer function  istato()
c-----------------------------------------------------------------------
c
c.... test if output file exists
c-----------------------------------------------------------------------
c
      USE comfil
      USE istat
      logical lout
      if(newf.ne.' ') then
        fout = newf
        newf = ' '
      endif
      wdo = 'new'
      inquire(file=fout,exist=lout)
      if(lout) wdo = 'exists'
      istato=1
      end
c
      integer function  istatr()
c.... test if restart(r) file exists
      USE comfil
      USE istat
      logical lres
      if(newf.ne.' ') then
        fres = newf
        newf = ' '
      endif
      wdr = 'new'
      inquire(file=fres,exist=lres)
      if(lres) wdr = 'exists'
      istatr=1
      end
c
      integer function  istatw()
c.... test if restart(w) file exists
      USE comfil
      USE istat
      logical lsav
      if(newf.ne.' ') then
        fsav = newf
        newf = ' '
      endif
      wds = 'new'
      inquire(file=fsav,exist=lsav)
      if(lsav) wds = 'exists'
      istatw=1
      end
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c.... FTN Graphik - Routinen
      subroutine vga
cww      call vga@
      end
c
      subroutine set_palette(ireg,ival)
      integer*2 ireg,ival
cww      call set_palette@(ireg,ival)
      end
c
      subroutine get_key(key)
      integer*2 key
cww      call get_key@(key)
      end
c
      subroutine get_mouse_position(ih2,iv2,bstat2)
      integer*2 ih2,iv2,bstat2
cww      call get_mouse_position@(ih2,iv2,bstat2)
      end
c
      subroutine set_graphics_selection(n1)
cww      include <windows.ins>
      integer*4 n1
cww      call set_graphics_selection@(n1)
      end
c
      subroutine create_polygon(xp,yp,n,handle,ierr)
      end
c
      subroutine fill_polygon(handle,icol,ierr)
      end
c
      subroutine delete_polygon_definition(handle,ierr)
      end
c
      subroutine set_text_attribute(font,size,rot,italic)
      integer*2 font
      real*4 size,italic,rot
cww      call set_text_attribute@(font,size,rot,italic)
      end
C
      integer*4 function iresize()
      iresize = 1
      return
      end
c
      subroutine pwmenu(i)
      end
c
      subroutine pwinhelp(name,wd,nwd,wrd)
      end
c
c----------------------------------------------------------------------
c.... dummy routinen UNIX
c
c----------------------------------------------------------------------
      subroutine getenv(i1,i2)
      end
c
c
c----------------------------------------------------------------------
c.... dummy routinen PHIGS,GKS
c
c----------------------------------------------------------------------
c.... dummy routinen PHIGS
c
      subroutine gpassw(i1,i2)
      end
c
      subroutine gplcmo(i1,i2,i3,i4)
      end
c
      subroutine gprqlc(i1,i2,i3,i4,i5)
      end
c
      subroutine gpqlc(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13)
      end
c
      subroutine gpinlc(i1,i2,i3,i4,i5,i6,i7,i8)
      end
c
      subroutine gptx2(i1,i2,i3)
      end
c
c----------------------------------------------------------------------
c.... dummy routinen GKS
c
      subroutine ginlc(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12)
      end
c
      subroutine gtx(i1,i2,i3)
      end
c
      subroutine grqlc(i1,i2,i3,i4,i5,i6)
      end
