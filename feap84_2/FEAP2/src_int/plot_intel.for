c-----------------------------------------------------------------------
c.... Additional programs for WINDOWS Intel
c
c-----------------------------------------------------------------------
c
      subroutine set_device(idev)
c-----------------------------------------------------------------------
c.... Purpose:  set type of compiler and graphcial device
c
c     Input    -
c
c     Output   idev,versn(1)
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE vdata
      integer idev

      versn(1) = 'VISTA/WIN7/WIN8 '
      idev=3

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine test_gen_op(ior,iow,finp,fout)
c-----------------------------------------------------------------------
c.... Purpose:  open the files for input and output  (from filnam!)
c
c     Input
c
c     Output
c
c     Comment: compilerdependent ->share
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------

      integer ior,iow,iorww
      character*229   finp,fout
      iorww=iabs(ior)
      open(unit=iorww,file=finp,status='unknown')
      open(unit=  iow,file=fout,status='unknown')
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine etimef(tt)
c-----------------------------------------------------------------------
c.... Purpose:  calculate time
c
c     Input    -
c
c     Output   tt time difference in sec.+1/100sec to start FEAP
c
c     Comment  tma: start time
c              tm: actual time
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE iofile
      USE timex
      real*4 tt,r,s,t,u
      integer ihrs,imins,isecs,ihunds

c...  time in seconds
      call GETTIM(ihrs,imins,isecs,ihunds)
      r = ihunds
      s = isecs
      t = imins
      u = ihrs
      tm = 3600.d0*u+60.d0*t+s+r*0.01
      if(jtime.eq.0) tma = tm
c.... delta t to begin of process
      tt = tm - tma
c.... for overnight calculations 24h=86400sec
      if (tt.lt.0) then
        tma = tma-86400
        tt  = tm -tma
      endif
      jtime=1
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine tstart
c-----------------------------------------------------------------------
c.... Purpose: full date/time info at start of execution
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE iofile

      integer*4 iday,imonth,iyear
      integer*4 ihrs,imins,isecs,ihunds

      CALL GETDAT(iyear,imonth,iday)
      CALL GETTIM(ihrs,imins,isecs,ihunds)

c.... plot time
                   write(iow,1000) iday,imonth,iyear,ihrs,imins,isecs
      if(ior.lt.0) write(*  ,1001) iday,imonth,iyear,ihrs,imins,isecs
1000  format(  1x,'FEAP-Calculation starts at  ',
     1        i2,'.',i2,'.',i4,6x,i2,':',i2,'.',i2)
1001  format(  9x,'FEAP-Calculation starts at  ',
     1        i2,'.',i2,'.',i4,6x,i2,':',i2,'.',i2)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine tend(lexst,md,mm,my)
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      logical lexst
      integer*4 id,im,iy,md,mm,my
      lexst=.true.

      CALL GETDAT(iy,im,id)

      if(iy.gt.my) then
        lexst=.false.
        return
      end if
      if(iy.lt.my) return
        if(im.gt.mm) then
          lexst=.false.
          return
        end if
        if(im.lt.mm) return
          if(id.gt.md) then
            lexst=.false.
            return
          end if
          if(id.le.md) return
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine startgr(isw)
c----------------------------------------------------------------------
c
c.... Purpose: open/close Graphic workstation INTEL
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      integer isw
#ifndef _BATCH_
      if(isw.eq.1) then

c....   Open Console Window
        call plotstwin

c....   Open Graphic Window
        call plotgrwin(1)

c....   some initial plot values
        call plotini

c....   Set Color tables
        call setcol(3,0)

      else if(isw.eq.2) then

c....   close Graphic Window

        call plotgrwin(2)

c....   close Console Window
        call clwclose(1,1,0)

      end if
#else
      call plotini
      call setcol(3,0)
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plopen
c-----------------------------------------------------------------------
c.... Purpose: open the plot device
c
c     Comment:
c     iclear = 0  open and clear window
c            = 1  open           window
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c
#ifndef _BATCH_
      use ifqwin
      USE iwinio
      USE pback
      USE pdata0
      USE pdata2
      USE pdata8
      USE pdatap
      USE pltran
      USE ppers
      implicit double precision (a-h,o-z)
      integer*2 icol
      integer*4 result
      result = SETACTIVEQQ(1)

      if(iclear.eq.0.and. .not. fopn) then
c
      end if
      result = SETBKCOLOR(16) ! white
      icol = 16
      if(iback.eq.1) then
         icol = 1
         result = SETBKCOLOR(32) ! black
      end if
      call pppcol(icol)
      if(iclear.eq.0) then
        call clear_screen_area
         istruc = 1
      end if
      ipan = 0
c...  set plotwindow = open
      fopn = .true.
c.... put up logo for feap  and border
      if(iclear .eq.0) then
        iclear = 1
        call plbord(1)
        call pfeap(0.98d0,0.02d0,0.240d0,2)
        if(kpers.eq.1.and.ibor.ge.0) call plview(qq,vmin,vmax)
        rots = sqrt(rotang(1)**2+rotang(2)**2+rotang(3)**2)
        if(rots.gt.0.0.and.kpers.eq.0.and.ibor.ge.0)
     +            call pltrot(tra,3,rotang,2)
      end if
#endif
      return
      end
c
      subroutine plclos
c-----------------------------------------------------------------------
c.... Purpose: close the plot device
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      use ifqwin
      USE pdata8
      USE pdatap
      implicit double precision (a-h,o-z)
      integer*4 result
c
      if(.not.fopn) return
      fopn = .false.
      result = SETACTIVEQQ(0) ! give active+control to textwindow
      result = FOCUSQQ(0)

c...  set new segment number
      istruc = istruc + 1
#endif
      return
      end
c

c----------------------------------------------------------------------
      subroutine clpan
c-----------------------------------------------------------------------
c.... Purpose: plot a filled polygon                                          |
c
c     Comments
c     poly = ixp,iyp = coordinates
c     ipan = No of points
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      use ifqwin
      USE pdata2
      USE pdatap
      USE pftn77
      USE plotter
      integer i
      integer*2 ipan2
      integer*4 ipan1
      integer*4 result
      real*8 xipan1,colps
      character wdummy*4
      type (xycoord) poly(400)

      ipan1=ipan
      do i = 1,ipan
        poly(i).xcoord = ixp(i)
        poly(i).ycoord = iyp(i)
      end do

c.... fill panel
      if(ipan.ge.3) then
        ipan2 = ipan
        result = POLYGON($GFILLINTERIOR,poly,ipan2)
        ipan = 0
      end if
c.... fill panel for PS
      if(nexte.gt.0.and.ihpgl.eq.1) then
        if(ipan1.ge.3) then
          xipan1 = ipan1
c.........find color
          colps = icps
          call psfile (3,colps,0.0d0,dummy,dummy,wdummy)
c.........plot panel
          call psfile (4,xipan1,1.0d0,dummy,dummy,wdummy)
c........ plot bounds of panel (not for disp etc.)
          if(ibps.eq.0) then
            call psfile (3, colps,1.0d0,dummy,dummy,wdummy)
            call psfile (4,xipan1,-1.d0,dummy,dummy,wdummy)
          end if
          ipan = 0
        end if
      end if
#endif
      return
      end
c
c-----------------------------------------------------------------------
      subroutine dplot(x,y,ipen,imove)
c-----------------------------------------------------------------------
c
c.... purpose: line drawing command
c
c     input:
c       x       - coordinate x1
c       y       - coordinate x2
c      ipen = 1 - begin panel at x,y
c      ipen = 2 - drawline    to x,y
c      ipen = 3 - move        to x,y
c      imove    - x=x+deltax and y=y+deltay
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE iwinio
      USE pdata2
      USE pdata11
      USE pdatap
      USE pftn77
      USE plotter
      implicit double precision(a-h,o-z)
      character wdummy*4
      integer idx,idy
      integer*2 ix12,iy12
      idx = iwxgs*0.75
      idy = iwygs
      ipan1 = ipan

c.... INTEL SCREEN
      jx1 =  x*idx
      jy1 =  y*idy
      if(imove.eq.1) then
        jx1 = (deltax+x)*idx
        jy1 = (deltay+y)*idy
      end if
      if (ipen .eq. 1) then
        ixa    = jx1
        iya = idy - jy1
        ixp(1) = ixa
        iyp(1) = iya
        ipan  = 1
      else if(ipen .eq. 2) then
        if(ipan.eq.0) then
          ix12 = jx1
          iy12 = idy - jy1
          call draw_line(ixa,iya,ix12,iy12,icc)
          ixa = ix12
          iya = iy12
        else if(ipan.ge.1) then
          ipan = ipan + 1
          iy12 = idy - jy1
          ixp(ipan) = jx1
          iyp(ipan) = iy12
          ixa  = jx1
          iya  = iy12
        end if
      else if(ipen .eq. 3) then
        iy12 = idy - jy1
        ixa    = jx1
        iya    = iy12
      end if
      if(nexte.gt.0) then
c.... HPGL (res = 7600) and PS
        xx1 =  x
        yy1 =  y
        if(imove.eq.1) then
          xx1 = deltax+x
          yy1 = deltay+y
        end if
        if (ipen .eq. 1) then
          xxp(1) = xx1
          yyp(1) = yy1
          ipan1 = 1
        else if(ipen .eq. 2) then
          if(ipan1.eq.0) then
            xxp(1) = xxp(2)
            yyp(1) = yyp(2)
            xxp(2) = xx1
            yyp(2) = yy1
            if(ihpgl.eq.1) call   psfile(2,xxp(1),yyp(1),xxp(2),yyp(2)
     +        ,wdummy)
            if(ihpgl.eq.2) call hpglfile(2,xxp(1),yyp(1),xxp(2),yyp(2))
          else if(ipan1.ge.1) then
            ipan1 = ipan1 + 1
            xxp(ipan1) = xx1
            yyp(ipan1) = yy1
          end if
        else if(ipen .eq. 3) then
          xxp(1) = xxp(2)
          yyp(1) = yyp(2)
          xxp(2) = xx1
          yyp(2) = yy1
        end if
      end if
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine drawtxt(ipl,x,y,icol,nchh,istr,str)
c-----------------------------------------------------------------------
c     purpose:
c     plot text on screen
c
c     input:
c        ipl  = 0 not on psfile,1 = on psfile
c        x,y  = position (FEAP)
c       icol = color (FEAP)
c       nchh = height of character
c       istr = length of text
c        str  = text
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
c      USE pdata2
      USE pdatap
      USE pftn77
      USE plotter
      implicit double precision(a-h,o-z)
      integer*2 istr
      character*(*) str
c.... set color, for cont+text on PC
      call pppcol(icol)
c.... set position
      call dplot(x,y,3,0)
c.... set character height
      call pltsiz(nchh)  ! ww to do
c.... draw text
      call draw_text(str,ixa,iya,icc)
c
      if(nexte.gt.0.and.ipl.eq.1) then
          call hptext (x,y,str,ihpgl)
      end if
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plline(iln)
c-----------------------------------------------------------------------
c     purpose: set line style                                                    |
c
c     input:
c     iln      line number
c
c     output:
c              line style
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      use ifqwin
      USE pdata2
      USE plinet
      USE plotter
      implicit double precision (a-h,o-z)
      character wdummy*4
      integer iln
      integer*2 ilph(8),istatus,lstyle
      data ilph /
     1 #FFFF , ! PS_SOLID
     2 #EEEE , ! PS_DASH
     3 #ECEC , ! PS_DASHDOT
     4 #ECCC , ! PS_DASHDOTDOT
     5 #AAAA , ! PS_DOT
     6 0 ,
     7 0 ,
     8 0 /
      lstyle=ilph(iln)
      CALL SETLINESTYLE(lstyle)
      ilinp=iln
c.... postscript output
      x1 = iln
      if(ihpgl.eq.1.and.nexte.gt.0)
     +   call psfile (7,x1,dummy,dummy,dummy,wdummy)
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltback(idev,iback,r,g,b)
c-----------------------------------------------------------------------
c     purpose: set background color
c
c     input:
c     iback   - 0: black on white
c     iback   - 1: white on black
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      use ifqwin
      real*4 r,g,b
      integer *4 idev,iback
      integer*4 status
      status = SETACTIVEQQ(1)
      if(iback.eq.1) then
        status = SETBKCOLOR(0)  ! black
        status = SETCOLOR(2)   ! bright white
      else
        status = SETBKCOLOR(16)  ! bright white
        status = SETCOLOR(1)     ! black
      end if
      CALL CLEARSCREEN ($GCLEARSCREEN)
      status = SETACTIVEQQ(0)
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pppcol(icol)
c-----------------------------------------------------------------------
c     Purpose: set color
c-----------------------------------------------------------------------
c.... Color map    FEAP
c       1  | black (white)
c       2  | red
c       3  | blue
c       4  | cyan
c       5  | green
c       6  | magenta
c       7  | yellow
c       8  | navy
c       9  | chocolate
c      10  | coral
c      11  | dark green
c      12  | dark red
c      13  | purple
c      14  | orange
c      15  | black (white)
c      16  | white (black)
c-----------------------------------------------------------------------
c     17-30: 14 colors for disp,stre etc. from red to blue
c-----------------------------------------------------------------------
c     white (black) | 31 | 0
c     black (black) | 32 | 0
c-----------------------------------------------------------------------
c     FEAP colors in brackets for iback = 1(black background)
c     FEAP colors             for iback = 0(white background)
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      USE colmap
      use ifqwin
      USE pback
      USE pdata2
      USE pftn77
      USE plinet
      USE plotter
      implicit double precision(a-h,o-z)
      integer*2 icol
      real*8 hpcol
      character wdummy*4
      integer*2 icolor
cww     save j
cww      data j/0/

      if(idev.eq.3) then
        if(icol.gt.0) then
          j = icol
          j = max(1,min(32,j))
          icc = icolmap(j)
          if(iback.eq.1) then
            if(j.eq. 1) icc=16
            if(j.eq.15) icc=16
            if(j.eq.16) icc=32
            if(j.eq.31) icc=32
          endif
c....     Line color
          icolor = SETCOLOR(icc)
        end if

c....   Interior fill style 2=solid 4 = hatch  1-4
cww        call gpis(2)
c....   Fill area color(solid)
c....   Text color
cww        call gptxci(j)
c....   Text hight
cww        call gpchh(.020)

c....   Line typ
c        is set in SR plline and in common plinet
      end if
c.... set values for HPGL and PS
      if(nexte.gt.0) then
c...    limit to colors 1 -- 32
        jhp   = icol
        jhp   = max(1,min(32,jhp))
        hpcol = jhp
        icps  = hpcol
        if(ihpgl.eq.1) then
c....     plot text in black
          call psfile (3,hpcol,1.0d0,dummy,dummy,wdummy)
        end if
        if(ihpgl.eq.2) call hpglfile (3,hpcol,dummy,dummy,dummy)
        end if
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine setcol(isw,imono)
c-----------------------------------------------------------------------
c     purpose: define colors and set color tables
c
c     input:
c     imono       not used
c     isw   - 3   define colors
c           - 1   not used
c           - 2   not used
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c                         Color tables
c   Graphic |  1    |  15   |  16   | 17-30 |  31   |  32
c-----------------------------------------------------------------------
c   color   | white | white | black | col   | white | black
c   mono    | white | white | white | grey  | black | black
c----------------------------------------------------------------------
#ifndef _BATCH_
      USE colmap
      use ifqwin
      integer*2 i
      integer*4 k,icolmap_i(32),result

c.... color table RGB HEX see http://de.wikipedia.org/wiki/Webfarben,
c                             but there: BRG!!

c.... colors table
      data  icolmap_i /
     1  #000000 ,   ! black
     2  #0045FF ,   ! red
     3  #FF0000 ,   ! blue
     4  #FFFF00 ,   ! cyan
     5  #00FF00 ,   ! green
     6  #FF00FF ,   ! magenta
     7  #00FFFF ,   ! yellow
     8  #800000 ,   ! navy
     9  #1E69D2 ,   ! chocolate
     +  #507FFF ,   ! coral
     1  #006400 ,   ! dark green
     2  #00008B ,   ! dark red
     3  #800080 ,   ! purple
     4  #008CFF ,   ! orange
c----------------------------------
     5  #000000 ,   ! black
     6  #FFFFFF ,   ! white
c----------------------------------
     7  #00008B ,   ! dark red
     8  #0000FF ,   ! red
     9  #0045FF ,   ! orange red
     +  #008CFF ,   ! dark orange
     1  #20A5DA ,   ! goldenrod
     2  #00FFFF ,   ! yellow
     3  #2FFFAD ,   ! green yellow
     4  #00FF00 ,   ! lime
     5  #7FFF00 ,   ! spring green
     6  #FFFF00 ,   ! cyan
     7  #FFBF00 ,   ! deep sky blue
     8  #FF901E ,   ! dodger blue
     9  #CD0000 ,   ! medium blue
     +  #8B0000 ,   ! dark blue
c----------------------------------
     1  #FFFFFF ,   ! brightwhite
     2  #000000 /   ! black
c----------------------------------


      if(isw.eq.1) then

c....   nothing to do

      else if(isw.eq.2) then

c....   nothing to do

      else if(isw.eq.3) then

c....   define colours 1-32
        do i = 1,32
          result = REMAPPALETTERGB(i,icolmap_i(i))
        end do

c       define standard sequence 1-32
        do i = 1,32
          icolmap(i) = i
        end do
        end if
#endif
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine getlcor1(xo)
c-----------------------------------------------------------------------
c
c.... Purpose: input screen coordinates  of plot line
c
c     Inputs:
c     xo = coordinates of line points (screen/feap)
c
c     Outputs
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
#ifndef _BATCH_
      use ifqwin
      use ifwin
      USE iwinio
      USE pdata11
      implicit double precision(a-h,o-z)
      integer*4 ixm,iym
      integer*4 mouseevent,keystate, result
      integer*4 cursor,newcursor
      real*8 xo(4)

      ixm = 0
      iym = 0
      mouseevent = MOUSE$LBUTTONDOWN

c...  set mouse cursor to cross hair
      newcursor = LoadCursor(0_HANDLE, IDC_CROSS)
      cursor    = SetMouseCursor(newcursor)

c.... move cursor to point 1, then press left button
      write(*,2000) 1
      result = WAITONMOUSEEVENT (mouseevent,keystate,ixm,iym)
      xo(1) =      ixm/(iwxgs*0.75)
      xo(2) = 1. - iym/(iwygs*1.00)

c.... move cursor to point 2 then press left button
      write(*,2001) 2
      result = WAITONMOUSEEVENT (mouseevent,keystate,ixm,iym)
      xo(3) =      ixm/(iwxgs*0.75)
      xo(4) = 1. - iym/(iwygs*1.00)

c.... modify coordinates do to move-option
      xo(1) = xo(1) - deltax
      xo(3) = xo(3) - deltax
      xo(2) = xo(2) - deltay
      xo(4) = xo(4) - deltay

c...  set mouse cursor to standard arrow
      newcursor = LoadCursor(0_HANDLE, IDC_ARROW)
      cursor    = SetMouseCursor(newcursor)

      return
2000  format(' mark ',i2,'. line point with left  mouse button')
2001  format(' mark ',i2,'. line point with left  mouse button')
#endif
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltzm(xpl)
c-----------------------------------------------------------------------
c
c.... Purpose: Graphical Input of bounds for ZOOM, see the clip macro
c
c     Output:
c       point 1 (x_1 = xpl(1), x_2 = xpl(2))
c       point 2 (x_1 = xpl(3), x_2 = xpl(4))
c
c     Comment:
c       acts only in 1 - 2 plane
c       valid also if move macro is acting on screen
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      use ifqwin
      use ifwin
      USE iwinio
      USE pdata11
      implicit double precision(a-h,o-z)
      integer*4 ixm,iym
      integer*4 mouseevent,keystate, result
      integer*4 cursor,newcursor
      dimension xpl(4)
c
      call plopen
      ixm = 0
      iym = 0
      mouseevent = MOUSE$LBUTTONDOWN

c...  set mouse cursor to cross hair
      newcursor = LoadCursor(0_HANDLE, IDC_CROSS)
      cursor    = SetMouseCursor(newcursor)

c.... get clip point 1
      write(*,2000) 1
      result = WAITONMOUSEEVENT (mouseevent,keystate,ixm,iym)
      x1 =      ixm/(iwxgs*0.75)
      x2 = 1.-  iym/(iwygs*1.00)
c.... including move
      xpl(1) = x1 - deltax
      xpl(2) = x2 - deltay
c.... control plot
      xp1 = x1
      xp2 = x2

c.... get clip point 2
      write(*,2000) 2
      result = WAITONMOUSEEVENT (mouseevent,keystate,ixm,iym)

      x1 =      ixm/(iwxgs*0.75)
      x2 = 1.-  iym/(iwygs*1.00)
c.... including move
      xpl(3) = x1 - deltax
      xpl(4) = x2 - deltay

c...  set mouse cursor to standard arrow
      newcursor = LoadCursor(0_HANDLE, IDC_ARROW)
      cursor    = SetMouseCursor(newcursor)

c.... plot actual rectangle (like rubberband, only control)
      call pppcol(7)
      call dplot(xp1 , xp2 ,3,0)
      call dplot(x1  , xp2 ,2,0)
      call dplot(x1  , x2  ,2,0)
      call dplot(xp1 , x2  ,2,0)
      call dplot(xp1 , xp2 ,2,0)

      return

2000  format(' mark ',i2,'. clip point with left  mouse button')
#endif
      end
c
c----------------------------------------------------------------------+
c
      subroutine pltpele(x,ix,ie,idis,nen,nen1,nie,ndm,numnp,numel)
c----------------------------------------------------------------------+
c
c     Purpose: Get number for an element using mouse
c     each left button click defines point,right button leads to end
c     including move,rot,iso,pers
c
c     Inputs:
c         x(ndm,*)   - Nodal coordinates
c        ix(nen1,*)  - Nodes on elements
c        idis        - local element numbers
c        nen         - no of Nodes   on elements(max)
c        nen1        - no of Nodes+3 on elements
c        numnp       - Number of nodal points
c        numel       - Number of elements
c
c     Outputs:
c        nz        - Selected element shown in text screen and plotted
c
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------+
c
#ifndef _BATCH_
      use ifqwin
      use ifwin
      USE iwinio
      USE pdata1
      USE pdata11
      USE ppers
      implicit double precision(a-h,o-z)
      integer*4 ixm,iym
      integer*4 mouseevent,keystate, result
      integer*4 cursor,newcursor
      character*5 yy
      dimension ix(nen1,*),ie(nie,*),idis(*),iplt(40),xx(3),xl(3,27)
      real*8 x(ndm,*),xg(3)
c
      dx1 = .005/scale
c
      nzo= 1
c.... set start position of cursor
      ixm = 0
      iym = 0

c...  set mouse cursor to cross hair
      newcursor = LoadCursor(0_HANDLE, IDC_CROSS)
      cursor    = SetMouseCursor(newcursor)

c.... set mouse event
      mouseevent = MOUSE$LBUTTONDOWN.or.MOUSE$RBUTTONDOWN


      write(*,*) 'LEFT Button: Get Element Numbers, RIGHT Button: End'

c.... first event
      result = WAITONMOUSEEVENT (mouseevent,keystate,ixm,iym)

c.... find point
10    x1 =     ixm/(iwxgs*0.75)
      y1 = 1.- iym/(iwygs*1.00)
c.... screen coordinates of point including move
      x1 = x1 - deltax
      y1 = y1 - deltay

c     Find closest center node of an element
      nz = 1
      xm = 1e6
      do n = 1,numel
        new = idis(n)
        if(iplma(ix(nen1,new)).eq.0)    goto 20 ! only if matn
        ma = ix(nen1,new)
        xx(1) = 0.0
        xx(2) = 0.0
        xx(3) = 0.0
        jj = 0
        call pltord(ie(nie-1,ix(nen1,n)),nn,iplt)
        nn = max(1,nn-1)
        do  i = 1,nn
          j  = iplt(i)
          if(j.gt.0 .and. j.le.nen) then
            ii = ix(j,n)
            if(ii.gt.0) then
              jj = jj + 1
              xx(1) = xx(1) + x(1,ii)
              xx(2) = xx(2) + x(2,ii)
              if(ndm.ge.3) xx(3) = xx(3) + x(3,ii)
            end if
          end if
        enddo
        if(jj.gt.0) then
c...    coordinates of center
          xx1 = xx(1)/jj
          xx2 = xx(2)/jj
          xx3 = xx(3)/jj
        end if
c....   perform perspective tranformation
        if(kpers.eq.1) then
          xg(1) = xx1
          xg(2) = xx2
          xg(3) = xx3
          call perspj(xg,xg,3,3,1,errv)
          xx1 = xg(1)
          xx2 = xg(2)
          xx3 = xg(3)
        end if
c....   screen coordinates of center point
        s1 = scale*(xx1 + xx1 -sx(1)) + s0(1)
        s2 = scale*(xx2 + xx2 -sx(2)) + s0(2)
c....   if isometric recompute values
        if(iso) then
          xmul = 2.*max(dx(1),dx(2))/(dx(1)+dx(2))
          te   = (s0(1) + 0.5*(s1 - s2))*xmul
c....     again screen coordinates of center point
          s2   = (0.2885*(s1 + s2))*xmul + scale*(xx3+xx3)
          s1   = te - 0.1
        end if
c....   distance (screen coordinates)
        ym = (s1-x1)**2+(s2-y1)**2
        if(ym.lt.xm) then
          xm = ym
          nz = n
        end if
20      continue
      end do
c.... print and plot element number (only if changes occur!!)
      if(nz.ne.nzo) then
c....   find center again and coordinates of all node
        call pzero(xl,81)
c        new = idis(nz)
c        ma = ix(nen1,new)
        ma = ix(nen1,nz)
        write(*,*) 'Element NO = ',nz,'  Material NO = ',ma
        xx(1) = 0.0
        xx(2) = 0.0
        xx(3) = 0.0
        jj = 0
        call pltord(ie(nie-1,ix(nen1,nz)),nn,iplt)
        nn = max(1,nn-1)
        do  i = 1,nn
          j  = iplt(i)
          if(j.gt.0 .and. j.le.nen) then
            ii = ix(j,nz)
            if(ii.gt.0) then
              jj = jj + 1
c
              xl(1,j) =  x(1,ii)
              xl(2,j) =  x(2,ii)
              if(ndm.eq.3) xl(3,j) = x(3,ii)
c
              xx(1) = xx(1) + x(1,ii)
              xx(2) = xx(2) + x(2,ii)
              if(ndm.eq.3) xx(3) = xx(3) + x(3,ii)
c
            end if
          end if
        enddo
c....   plot element
        call pppcol(5)
        call plot9(ie(nie-1,ma),ix(1,nz),xl,3,nen,2)
        if(jj.gt.0) then
c....   coordinates of center
          xx1 = xx(1)/jj
          xx2 = xx(2)/jj
          xx3 = xx(3)/jj
        end if
        xx1 = xx1 - dx1*5
        xx2 = xx2 - dx1
c....   plot node number
        call plotl(xx1,xx2,xx3,3)
        call pppcol(3)
        call plabl(nz)
c
        nzo = nz
      end if

      result = WAITONMOUSEEVENT (mouseevent,keystate,ixm,iym)

      if(result.eq.1) go to 10

c...  set mouse cursor to standard arrow
      newcursor = LoadCursor(0_HANDLE, IDC_ARROW)
      cursor    = SetMouseCursor(newcursor)
#endif
      return
      end
c
c----------------------------------------------------------------------+
c
      subroutine pltpnod(x0,x,b,ds,anggl,ix,nen1,ndm,ndf,numnp,flstre,
     1 k2,k3,ipola)
c----------------------------------------------------------------------+
c
c     Purpose: Get number for a node using mouse
c     and plot available values for disp,stre,reac
c     each left button click defines point,right button leads to end
c     including move,rot,iso,pers
c
c     Inputs:
c        x0(ndm,*) - original coordinates      |
c        x(ndm,*)  - Nodal coordinates in def. state x=x+c*b with rot
c        b(ndf,*)  - Nodal displacements
c       ds(numnp,*)- Nodal stresses
c      anggl(numnp)- Nodal angles
c       ix(nen1,*) - Nodes on Elements
c       k2         - .ne.0 additional results
c       k3         - .ne.0 write on file
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------+
c
#ifndef _BATCH_
      use ifqwin
      use ifwin
      USE iofile
      USE iwinio
      USE pdata1
      USE pdata11
      USE ppers
      USE strnam
      implicit double precision(a-h,o-z)
      logical flstre
      dimension ix(*)
      integer*4 ixm,iym
      integer*4 mouseevent,keystate, result
      integer*4 cursor,newcursor
      real*8 x0(ndm,*),x(ndm,*),b(ndf,*),ds(numnp,*),anggl(numnp),
     1       xg(3),bl(7)
c
      ista = iabs(istv)
c
      nzo= 1
c.... set start position of cursor at midpoint
      ixm = 0
      iym = 0

c.... set mouse event
      mouseevent = MOUSE$LBUTTONDOWN.or.MOUSE$RBUTTONDOWN

c.... set mouse cursor to cross hair
      newcursor = LoadCursor(0_HANDLE, IDC_CROSS)
      cursor    = SetMouseCursor(newcursor)

      write(*,*) 'LEFT Button: Get Node Numbers, RIGHT Button: End'


c.... first event
      result = WAITONMOUSEEVENT (mouseevent,keystate,ixm,iym)

c.... find point
10    x1 =     ixm/(iwxgs*0.75)
      y1 = 1.- iym/(iwygs*1.00)
c.... screen coordinates of point including move
      x1 = x1 - deltax
      y1 = y1 - deltay
c     Find closest node
      nz = 1
      xm = 1e6
      do n = 1,numnp
        if(iplmano(ix,n,nen1).eq.0)  goto 100  ! only if matn
        xx1 = x(1,n)
        xx2 = x(2,n)
        if(ndm.eq.3) xx3 = x(3,n)
c....   perform perspective tranformation
        if(kpers.eq.1) then
          xg(1) = xx1
          xg(2) = xx2
          xg(3) = xx3
          call perspj(xg,xg,3,3,1,errv)
          xx1 = xg(1)
          xx2 = xg(2)
          xx3 = xg(3)
        end if
c....   screen coordinates of point
        s1 = scale*(xx1 + xx1 -sx(1)) + s0(1)
        s2 = scale*(xx2 + xx2 -sx(2)) + s0(2)
c....   if isometric recompute values
        if(iso) then
          xmul = 2.*max(dx(1),dx(2))/(dx(1)+dx(2))
          te   = (s0(1) + 0.5*(s1 - s2))*xmul
c....     again screen coordinates of point
          s2   = (0.2885*(s1 + s2))*xmul + scale*(xx3+xx3)
          s1   = te - 0.1
        end if
c....   distance (screen coordinates)
        ym = (s1-x1)**2+(s2-y1)**2
        if(ym.lt.xm) then
          xm = ym
          nz = n
        end if
100     continue
      end do
c.... print and plot point (only if changes occur!!)
      if(nz.ne.nzo) then
                    write(  *,*) 'NODE = ',nz
        if(k3.ne.0) write(iow,*) 'NODE = ',nz
        call pltnod(x,ix,nen1,ndm,numnp,-nz,7)
c....   print available values
        if(k2.ne.0) then
c....     coordinates
                      write(  *,2000) (i,' coord',i=1,ndm)
          if(k3.ne.0) write(iow,2000) (i,' coord',i=1,ndm)
          do i = 1,ndm
            xg(i) = x0(i,nz)
          end do
                      write(  *,2001) (xg(i),i=1,ndm)
          if(k3.ne.0) write(iow,2001) (xg(i),i=1,ndm)
c....     displacements
                     write(  *,2002) (i,' displ',i=1,ndf)
          if(k3.ne.0)write(iow,2002) (i,' displ',i=1,ndf)
c....     global displacements
          ndff=min(ndf,7)
          do i = 1,ndff
            bl(i) = b(i,nz)
          end do
          if(ipola.ne.0) then
c....       in case of polar coordinates
            angl = anggl(nz)
            call panglb1(xg,bl,angl,ndm,ndff,ipola)
          end if
                      write(  *,2001) (bl(l),l=1,ndff)
          if(k3.ne.0) write(iow,2001) (bl(l),l=1,ndff)
c....     stresses
          if(flstre) then
            if(strsus(1).eq.'               ') then
                          write(  *,2003)(i,i=1,ista)
              if(k3.ne.0) write(iow,2003)(i,i=1,ista)
            else
                          write(  *,2004)(strsus(i),i=1,ista)
              if(k3.ne.0) write(iow,2004)(strsus(i),i=1,ista)
            end if
                        write(  *,2005) (ds(nz,l),l=1,ista)
            if(k3.ne.0) write(iow,2005) (ds(nz,l),l=1,ista)
          end if
        end if
        nzo = nz
      end if

      result = WAITONMOUSEEVENT (mouseevent,keystate,ixm,iym)

      if(result.eq.1) go to 10

c...  set mouse cursor to standard arrow
      newcursor = LoadCursor(0_HANDLE, IDC_ARROW)
      cursor    = SetMouseCursor(newcursor)

      return

2000  format('  Coordinates ',/6(i6,a6):/(6(i6,a6)))
2001  format(1p6e12.5:/(1p6e12.5))
2002  format('  Displacements',/6(i6,a6):/(6(i6,a6)))
2003  format('  Stresses'/5(i7,'-value')/(5x,5(i7,'-value')))
2004  format('  Stresses'/5(1x,a15)/(5x,5(1x,a15)))
2005  format(1p5e16.5/(1p5e16.5))
#endif
      end
c
c-----------------------------------------------------------------------
c
      subroutine doargs1(nargs)
c-----------------------------------------------------------------------
c
c.... Purpose: count number of arguments on command line
c
c     Output:  nargs - number of arguments
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      integer nargs
      nargs = command_argument_count()
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine doargs2(argv,nchars,i)
c-----------------------------------------------------------------------
c
c.... Purpose: read one argument on command line
c
c     Output:
c       i      - position
c       nchars - not used, should be no. of characters
c       argv   - argument
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      integer i
      character*229 argv
      call getarg(i,argv)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine getparaf(fpath,fini)
c-----------------------------------------------------------------------
c.... Purpose: get fcfg and ini-name
c
c     Output:
c       fpath - path to feapwini.ini
c       fini  - feapwin.ini
c
c     W. Wagner IBS KIT 01/14
c-----------------------------------------------------------------------
      use ifport
      integer*4 result
      character*229 fpath
      character*11 fini
      result = GETENVQQ ('FCFG',fpath)
      fini = 'feapwin.ini'
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine getparap(ppath)
c-----------------------------------------------------------------------
c.... Purpose: get pcfg
c
c     Output:
c       fpath - path to feap directory
c
c     W. Wagner IBS KIT 01/14
c-----------------------------------------------------------------------
      use ifport
      integer*4 result
      character*229 ppath
      result = GETENVQQ ('PCFG',ppath)
      return
      end
c
c-----------------------------------------------------------------------
      subroutine plmousec(istat)
c-----------------------------------------------------------------------
c
c.... Purpose: get state of mouse
c
c     Output:
c      istat  - button status: 0=no button,1=left button, 2=right button
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
#ifndef _BATCH_
      use ifqwin
      integer*4 istat
      integer*4 ixm,iym
      integer*4 mouseevent,keystate, result

      mouseevent = MOUSE$LBUTTONDOWN.or.MOUSE$RBUTTONDOWN
     +            .or.MOUSE$LBUTTONDBLCLK

      result = WAITONMOUSEEVENT (mouseevent,keystate,ixm,iym)

      if(result.eq.1)istat = 1
      if(result.eq.8)istat = 2
      if(result.eq.4)istat = 3 ! double click: result_1=1, result_2=4
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltsiz(nn)
c-----------------------------------------------------------------------
c
c.... Purpose: set character height of text
c
c     Input:
c       nn    - Number of size
c
c     Output:
c             - isize and xx
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE hptext1
      implicit double precision (a-h,o-z)
      real*4 xx
      xx = 0.15 * nn + 0.85
      call set_text_attribute(1,xx,0,0)
      isize = nn
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine start_process(fcis)
c-----------------------------------------------------------------------
c
c.... Purpose: start an external process and wait until return
c
c     Inputs:
c       fcis  - name of program with path and options
c
c     Outputs
c       i     - return value
c               0 or >0 successful
c               -1 indicates an error condition
c
c     W. Wagner BS KIT 11/13
c
c-----------------------------------------------------------------------
      use ifport
      integer*4 i,errnum
      character*229 fcis
c
      i = system(fcis)

      if (i .eq. -1) then
        errnum = ierrno( )
       write(*,*) 'error in start_process with',fcis,
     +            'error number', errnum
      end if

      return
      end
c
c----------------------------------------------------------------------
c
      subroutine start_pprocess(fcis,fcis1)
c-----------------------------------------------------------------------
c
c.... Purpose: start an external process and return immediately
c              ### no such routine found for INTEL
c              ### used: runqq has 2 parameter but waits!!
c
c     Inputs:
c       fcis  - name of program
c       fcis1 - additional data
c
c     Outputs:
c       ifail - control parameter
c
c     Comments:
c     1) feapwin.ini:
c     fadobe Acrobat Reader (ohne ""!)
c     C:\Program Files (x86)\Adobe\Acrobat 10.0\Acrobat\Acrobat.exe
c     oder ...\Adobe\Acrobat 10.0\Acrobat\AcroRd32.exe
c     oder ...\Foxit Reader\Foxit Reader.exe
c
c     2) FEAP geht nicht weiter, sondern
c
c     3) PDF muss geschlossen werden!
c
c     4) GOTOACRO geht derzeit nicht!
c
c     W. Wagner BS KIT 12/13
c-----------------------------------------------------------------------
      use ifport
      integer*2 result
      integer*4 errnum
      character fcis*229,fcis1*229

c
      write(*,*) 'Please close PDF to return to FEAP'

      result = runqq(fcis,fcis1)

      if (result .eq. -1) then
        errnum = ierrno( )
       write(*,*) 'error in start_process with',fcis,fcis1,
     +            'error number', errnum
      end if

      return
      end
c
c
      subroutine start_pprocess_adobe(fcis,fcis1)
c-----------------------------------------------------------------------
c
c.... Purpose: start an external process and return immediately
c
c     Inputs:
c       fcis  - name of program
c       fcis1 - additional data
c
c     Outputs:
c       i     - return value
c               0 successful
c              -1 indicates an error condition
c
c
c     Comments: Starts Adobe Reader correct but only once! Adobe problem!
c
c     at present not used in FEAP
c
c     W. Wagner IBS KIT 11/12
c     D. Heller     TUD
c-----------------------------------------------------------------------
#ifndef _LINUX_
      use Kernel32
      implicit none
      integer(kind=2) i,ipos
      integer(bool) :: res
      character(len=229) :: fcis, fcis1

      type( t_process_information ) tPI
      type( t_startupinfo ) tSI

      call zeromemory(loc(tSI), sizeof(tSI))
      tSI%cb = sizeof(tSI)

c.....change \ in \\
      call enrich_fcis(fcis,229)
      call enrich_fcis(fcis1,229)

c ....add closing C
      i = ipos(fcis,229)
      fcis(i+1:i+1) = char(0)
      i = ipos(fcis1,229)
      fcis1(i+1:i+1) = char(0)

c.... start process
      res = CreateProcess(fcis,fcis1,0,0,1_BOOL,0_DWORD,0,0,tSI,tPI)

      if (res .eq. 0) then
        write(16,*) "Error! Could not create process.",GetLastError()
      end if

      res = CloseHandle(tPI%hProcess)
      res = CloseHandle(tPI%hThread)
#endif
      return
      end

c-----------------------------------------------------------------------

      subroutine enrich_fcis(fcis,nn)
c----------------------------------------------------------------------
c
c      Purpose: enrich in fcis string \ to \\ for use in C
c
c      Inputs:
c         fcis   with \
c
c      Outputs:
c         fcis   with \\
c
c      Comments: FTN    not necessary, but do not disturb
c                INTEL  necessary
c
c     W. Wagner BS KIT 11/12
c----------------------------------------------------------------------
      implicit none
      integer i,k,nn,n1,ipos
      character*1 fcis(nn),fcis1(nn)

      n1 = ipos(fcis,nn)

      fcis1=' '
      k=1
      do i=1,n1
        if(fcis(i).eq.'\') then
          fcis1(k)  ='\'
          fcis1(k+1)='\'
          k=k+2
        else
          fcis1(k)=fcis(i)
          k=k+1
        end if
      if(k.gt.nn) then
        stop 'fcis1 too long! in enrich_fcis'
      end if
      end do

      fcis=' '
      do i=1,k
        fcis(i)=fcis1(i)
      end do

      return
      end

c----------------------------------------------------------------------
c
      subroutine start_process_fe2(fcis,fcis1)
c-----------------------------------------------------------------------
c
c.... Purpose: start FEAP process inside from FEAP
c
c     Inputs:
c       fcis  - FEAP
c       fcis1 - additional data
c
c     Outputs
c       ifail - control parameter
c
c     Comment
c
c     W. Wagner BS KIT 10/10
c-----------------------------------------------------------------------
      use ifport
      integer*2 result
      integer*4 errnum
      character fcis*229,fcis1*600
c
      result = runqq(fcis,fcis1)

      if (result .eq. -1) then
        errnum = ierrno( )
       write(*,*) 'error in start_process with',fcis,fcis1,
     +            'error number', errnum
      end if

      return
      end
c
      subroutine get_number_proc(nproc)
c----------------------------------------------------------------------
c
c.... Purpose: get number of processors                            
c
c     Input:  
c       nproc - number of processors to be used
c
c     Output:  
c       mproc - number of processors available 
c       nproc - number of processors to be used 
c
c     W. Wagner BS UKA 09/14
c----------------------------------------------------------------------
      USE omp_lib
      integer*4 nproc,mproc

      mproc = omp_get_max_threads()  

      if(nproc.eq.0) then
        nproc=mproc              
      else
        nproc=min(nproc,mproc)    
        call omp_set_num_threads(nproc)
      end if

      return
      end
c



















c################################# to do


c-----------------------------------------------------------------------
c
c
      subroutine prin_inp(iname,filein)
c----------------------------------------------------------------------
c.... print existing input files, routine not used
c----------------------------------------------------------------------
      character iname*229,filein*229
      nc   = ipos(filein,229)
        write(*,2000) filein
      return
2000  format(/'ERROR - Specified input file ',a,/,
     +        ' does not exist, reinput name.'/)
      end
c
