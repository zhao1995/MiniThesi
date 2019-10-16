      subroutine hpglfile (ifxn,x1,y1,x2,y2)
c-----------------------------------------------------------------------
c      Purpose: writes a hpgl output file
c
c      Input:
c        ifxn   - 1 : open the plot file x1=filenumber fpgl=filename
c                 2 : draw line from (x1,y1) to (x2,y2)
c                 3 : select the line color   x1
c                 4 : dummy
c                 5 : close the plot file
c       x1      - set with respect to ifxn
c       y1      -
c       x2      -
c       y2      -
c
c      Output:
c       data on PS-File
c
c      Comment:
c      the dimension of x and y is:
c      (0,0) lower left  corner, (1,1) upper right corner
c      knebel oct. 94
c-----------------------------------------------------------------------
      USE hpgl1
      USE lplot1
      USE plotter
      implicit double precision(a-h,o-z)
      character*12  string1,string2,y*1
      integer*2 icolhp(8)
      logical lflg
      data lflg/.false./

c
c.... "xgupcm" is x graphics units per centimeter (set for din a4)
c
      data xgupcm,ygupcm /7600.0,7600.00/
c.... color-palette for feap->vga->hpgl viewer
      data icolhp /1,4,7,3,8,6,2,5/
c
      if(iprin.eq.0) return

      go to (100,200,300,400,500) ifxn
c
c.... open the  device
100   continue
      lun = x1
      inquire(file=fpgl,exist=lflg)
      if(lflg) then
        write(*,2000) fpgl
        read (*,1000) y
        if(y.eq.'y' .or. y.eq.'Y') go to 101
        return
101   continue
      end if
      open (unit=lun,file=fpgl,status='unknown')
      write(lun,'(a)') 'IN;SP1;'
      lhpgl = .true.
      return
c
c.... draw line from  x1,y1  to  x2,y2
200   continue
      ixposn = xgupcm*x1+0.5
      iyposn = ygupcm*y1+0.5
      call gdhpglconvert ( ixposn,iyposn,string1 )
      ixposn = xgupcm*x2+0.5
      iyposn = ygupcm*y2+0.5
      call gdhpglconvert ( ixposn,iyposn,string2 )
      write(lun,'(a)') 'PU;PA'//string1//'PD;PA'//string2
      return
c
c.... select current drawing color
300   continue
      icolo = x1
      if (icolo .le. 0 .or. icolo .gt. 8) then
        icolor = 1
      else
       icolor= icolhp(icolo)
      end if
      string1=' '
      call nnumstr(icolor,string1,idummy)
      write(lun,'(a)') 'PU;SP'//string1//';'
      return
c
400   return
c
c.... close the device
500   continue
      write(lun,'(a)') 'PU;SP0;'
      close (unit=lun)
      lhpgl = .false.
      return
c.... formats
1000  format(a1)
2000  format(' ** WARNING ** File ',a229,/,
     +       ' exists. Erase? (y or n) >',$)
c
      end
c
      subroutine gdhpglconvert(ix,iy,string)
c-----------------------------------------------------------------------
c
c      Purpose: converts the (x,y) pair into the proper hpgl string
c
c      Input:
c       ix       - coordinate x
c       iy       - coordinate y
c
c      Output:
c       string   - string of hpgl
c
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      character*(*) string
      string=' '
      call nnumstr(ix,string,iend1)
c     iend = leng(string)
      iend = iend1
      string(iend:iend) = ','
      call nnumstr(iy,string(iend+1:12),iend1)
      iend = iend+iend1
c     iend = leng(string)
      string(iend:iend) = ';'
      return
      end
c
      subroutine nnumstr(jval,bstrng,j)
c------------------------------------------------------------------------
c
c      Purpose: converts "jval" to a string with no leading
c               or trailing spaces.
c
c      Input:
c       jval
c
c      Output:
c       bstring
c
c------------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      character*(*) bstrng
      character*10  number, num2
      character*1   bnum(10)
      equivalence (num2,bnum(1))
c
      write(number,11) jval
11    format(i9)
      read(number,12) num2
12    format(a10)
      j=1
      do 100 i=1,9
         if (bnum(i) .eq. ' ') go to 100
         bstrng(j:j)=bnum(i)
           j=j+1
100   continue
      return
      end
c
      subroutine hptext(xx,yy,bstrng,ihpgl)
c---------------------------------------------------------------
c
c      Purpose: does the hpgl and postscript text plotting
c
c      Input:
c       xx     - x-position of the first character
c       yy     - y-position of the first character
C       ihpgl  - 1 => postscript
C              - 2 => hpgl
c       bstrng - text to plot
c
c       Output:
c       HPGl/PS-data
c
c      knebel oct. 94
c---------------------------------------------------------------
      USE hptext1
      implicit double precision(a-h,o-z)
      character*(*) bstrng,pstrng*80
      logical  lmove
      integer*2 bxy,index
      common /gccpar/  csize, ccos, csin, xoff, yoff
      common /gcfont/  indx(96),bxy(1556)
      common /psplot1/ scl,xxoff,yyoff
      if(ihpgl .eq.2 ) then
c.... hpgl-text
c
      ccos = 1.d0
      csin = 0.d0
c
      xvpos = xx
      yvpos = yy
c.... csize is the  size of the characters
c.... isize is input via plot-macro [text,,isize]  default = 1
      csize = 0.00175*isize
c
      nbyte = 0
      ii = len(bstrng)
      do 100 i = 1,ii
      nbyte = nbyte + 1
c
      xoff = xvpos
      yoff = yvpos
c.... get the character to stroke
      iichar = ichar(bstrng(nbyte:nbyte))
      if (iichar .eq. 0) return
c
c.... stroke the character

c.... space fill all non-printing and non-defined with space size of
c.... a capital "A".
c
      jchar = (iichar-32)
c
c.... stroke this character
c
      if (jchar .eq. 0) then
         index = 0
      else
         index = indx(jchar)
      end if
      idone = indx(jchar+1)
c
c.... first position is an assumed move
c
      lmove = .true.
c
c.... get the scaled and rotated next position on the character
c
110   continue
      if (bxy(index).eq.-64) then
         lmove = .true.
         index = index + 1
         go to 110
      end if
c
      x=bxy(index)
      y=bxy(index+1)
c
c     the following does the character sizing and rotation.
c
      xs1 = x*csize
      ys1 = y*csize
      dx = ccos*xs1 + csin*ys1 + xoff
      dy = ccos*ys1 - csin*xs1 + yoff
      index = index + 2
      xvpos = dx
      yvpos = dy
        if (lmove) then
           x2 = dx
           y2 = dy
        else
           x1 = x2
           y1 = y2
           x2 = dx
           y2 = dy
           if(ihpgl.eq.2) call hpglfile ( 2,x1,y1,x2,y2 )
c.kneb     if(ihpgl.eq.1) call   psfile ( 2,x1,y1,x2,y2,dummy )
        end if
        lmove = .false.
        if (index .lt. idone) go to 110
c
c...  all done with the character, move to next character position
c
      xs1 = 9.0*csize
      ys1 = 0.0
      dx = ccos*xs1 + csin*ys1 + xoff
      dy = ccos*ys1 - csin*xs1 + yoff
      xvpos = dx
      yvpos = dy
      x2 = dx
      y2 = dy
100   continue
      else if(ihpgl.eq.1) then
c....   postscript-text
C       move-to
        xpp = scl*xx+xxoff + 10.
        ypp = scl*yy+yyoff - 10.
        pstrng = ' '
        write(pstrng,'(f5.1,1x,f5.1,1x,a2)')  xpp,ypp,'m '
        call gspst(pstrng,14)
        il=len(bstrng)
        ie=il
        do  i=il,1,-1
           if (bstrng(i:i) .ne. ' ') go to 101
           ie=ie-1
        end do
101     ib=0
        do  i=1,il
           if (bstrng(i:i) .ne. ' ') go to 102
           ib=ib+1
        end do
102     il = ie-ib
        pstrng = ' '
        pstrng(1:il+5)  = '('//bstrng(ib+1:ie)//') s '
        call gspst(pstrng,il+5)
      else
      end if
      return
      end
c
C
      subroutine psfile (ifxn,x1,y1,x2,y2,macro)
c-----------------------------------------------------------------------
c
c      Purpose: writes a postscript plot file
c
c      Input:
c       ifxn   - 1 : open the plot file x1=filenumber fpgl=filename
c                2 : draw line from (x1,y1) to (x2,y2)
c                3 : select the line color  x1   (y1 for mono)
c                4 : plot polygon x1=number of points
c                    points are in  xxp,yyp (common /plotter/)
c                    filled for y1 > 0 not filled for  y2 < 0
c                5 : close the plot file
c                6 : define colors and grey shading
c                7 : set line-style x1
c                8 : scale text by x1
c                9 : write feap macro
c       x1    -  with respect to ifxn
c       y1    -
c       x2    -
c       y2    -
c       macro - Name of FEAP plot macro
c
c      Output:
c
c
c      Comments:
c        the dimension of x and y is (0,0) lower left  corner
c                                    (1,1) upper right corner
c        colors:   1-16 = black(lines), grey(fill)  1-16 = color (lines)
c                 17-30 = grey (fill)              17-30 = color (fill)
c                 31-33 = black                    31-32 = black
c
c
c     knebel oct. 94
c-----------------------------------------------------------------------
      USE hpgl1
      USE lplot1
      USE pdata2
      USE plotter
      implicit double precision(a-h,o-z)
      character*25 coord
      character*4 macro,y*1
      logical lflg
      common /colini/ icolold
      common /ipsbuf/  ibuff,ibuf(481)
      common /psplot1/ scl,xoff,yoff
      data lflg/.false./
c
c     scl  : size of the picture
c     xoff : horizontal offset (movement of the entire picture)
c     yoff : vertikal     offset
c
      if(iprin.eq.0) return
      go to (100,200,300,400,500,600,700,800,900) ifxn
c
c.... open the device
100   continue
      lun = x1
      inquire(file=fpgl,exist=lflg)
      if(lflg) then
        write(*,2000) fpgl
        read (*,1000) y
        if(y.eq.'y' .or. y.eq.'Y') go to 101
        return
101   continue
      end if
      open (unit=lun,file=fpgl,status='unknown')
c
      write(lun,1000) '%!PS-Adobe-2.0'
      if(ilsc.eq.2) then
        write(lun,1000) '%%BoundingBox: 50 80 560 740'   !landscape
      else
        write(lun,1000) '%%BoundingBox: 40 40 570 450'   !portrait
      end if
      write(lun,1000) '%%EndComments'
      write(lun,1000) '1 setlinecap 1 setlinejoin'
      write(lun,1000) '0.4 setlinewidth   % Liniendicke'
      write(lun,1000) '/m {moveto} def '
      write(lun,1000) '/l {lineto} def '
      write(lun,1000) '/v {stroke newpath} def'
      write(lun,1000) '/f {fill newpath} def'
      write(lun,1000) '/s {show} def'
      write(lun,1000)
     +        '/sf {/Times-Roman findfont exch scalefont setfont} def'
      if(ilsc.eq.2) then
c.... landscape
       write(lun,1000) '%for portrait delete the following 3 lines'
       write(lun,1000) '0 792  translate'
       write(lun,1000) '-90 rotate'
       write(lun,1000) '1.25 1.25 scale'
      else
c.... portrait
       write(lun,1000) '%for landscape uncomment the following 3 lines'
       write(lun,1000) '%0 792  translate'
       write(lun,1000) '%-90 rotate'
       write(lun,1000) '%1.25 1.25 scale'
       write(lun,1000) '%for seascape uncomment the following 3 lines'
       write(lun,1000) '%560 0  translate'
       write(lun,1000) '%90 rotate'
       write(lun,1000) '%1.25 1.25 scale'
      end if
cww   write(lun,1000) '8 sf    %Schriftgroesse auf Standardwert'
      write(lun,1000) '12 sf    %Schriftgroesse auf Standardwert'
c
      ibuff=0
      lps = .true.
c.... formats
1000  format(a)
2000  format(' ** WARNING ** File ',a229,/,
     +       ' exists. (PRIN,100 shows all plotfiles in use!)',/,
     +       ' Erase? (y or n) >',$)
      return
c
c.... draw line from  x1,y1  to  x2,y2
200   continue
c
c.... move to
c
      xpp = scl*x1+xoff
      ypp = scl*y1+yoff
      write(coord,'(f5.1,1x,f5.1,1x,a2)')  xpp,ypp,'m '
      call gspst(coord,14)
c
c.... draw to
c
      xpp = scl*x2+xoff
      ypp = scl*y2+yoff
      write(coord,'(f5.1,1x,f5.1,1x,a2)')  xpp,ypp,'l '
      call gspst(coord,14)
c
      return
c
c.... select current drawing color
300   continue
      icol = x1
      ityp = y1
      if(imono.ne.0.and.ityp.eq.1) icol = 33
c.... write only in case of new color (ityp not considered)
      if(icol.eq. icolold) return
      icolold=icol
      if(imono.ne.0) then
c       grey colors
cww        if(ityp.eq.1) icol = 33
        if(icol.le.9) write(coord,'(a3,i1,1x)') '  g',icol
        if(icol.gt.9) write(coord,'(a2,i2,1x)')  ' g',icol
      else if(imono.eq.0) then
c       colors
        if(icol.le.9) write(coord,'(a3,i1,1x)') '  c',icol
        if(icol.gt.9) write(coord,'(a2,i2,1x)')  ' c',icol
      end if
      call gspst(coord,5)
      return
c
c.... plot filled polygon
400   continue
      npts = x1
      nfill= y1
      call psfill (xxp,yyp,npts,nfill)
      return
c
c.... close the device
500   continue
      call clsbuf
      close (unit=lun)
      lps = .false.
      return
c
c.... define colors
c     #  1-16 = grey for standard  ,17-30 grey shading for fill 30-33 black
c        1-32 = black for lines
c     #  1-16 = standard colors,17-30 colors for fill
c
600   continue
c.... set values grey
      do i = 1,16
cww        clst(i)    = 0.15 + i*0.050d0
cww        clst(i+16) = 0.28 + i*0.045d0
        clst(i)    = 0.04 + (i-1)*0.07d0
        clst(i+16) = 0.04 + (i-1)*0.07d0
      end do
c.... set value black
      clst(31) = 0.0d0
      clst(32) = 0.0d0
      clst(33) = 0.0d0
      write(lun,1000)
     +'%define grey shading: 1-16=standard,17-30=fill,31-33=black'
      do i=1,33
        col  = clst(i)
        if(i.le.9) write(lun,6001) i,col
        if(i.gt.9) write(lun,6002) i,col
      end do
6001  format('/g',i1,' {  v ',f4.2,1x,'setgray } def')
6002  format('/g',i2, ' { v ',f4.2,1x,'setgray } def')
c.... set values color
      write(lun,1000)
     +'%define colors: 1-16=standard,17-30=fill,31-32=black'
      do i =1,32
        colr = colps(1,i)
        colg = colps(2,i)
        colb = colps(3,i)
        if(i.le.9) write(lun,6003) i,colr,colg,colb
        if(i.gt.9) write(lun,6004) i,colr,colg,colb
      end do
6003  format('/c',i1,' {  v ',3(f4.2,1x),'setrgbcolor } def')
6004  format('/c',i2, ' { v ',3(f4.2,1x),'setrgbcolor } def')
      return
c
c...  set line-style
700   continue
      iltyp = x1
      coord=' v '
      call gspst(coord,3)
      if(iltyp.eq.1) coord='[] 0            setdash ' ! solid
      if(iltyp.eq.2) coord='[1 2] 0         setdash ' ! dotted
      if(iltyp.eq.3) coord='[4 2 1 4] 0     setdash ' ! dash-dot
      if(iltyp.eq.4) coord='[4 2] 0         setdash ' ! short dash
      if(iltyp.eq.5) coord='[8 4] 0         setdash ' ! long dash
      if(iltyp.eq.6) coord='[1 1 1 1 4 6] 0 setdash ' ! dot-dot-dash
      if(iltyp.eq.7) coord='[4 2 8 2] 0     setdash ' ! short dash-long dash
      if(iltyp.eq.8) coord='[4 6] 0         setdash ' ! wide dash
      call gspst(coord,24)
      return
c
c...  scale text
800   continue
      if(x1.le.0.d0) return
      write(coord,'(a2,f5.1,a4)') ' v', 8.0*x1,' sf '
      call gspst(coord,11)
      return
c...  write feap macro
900   continue
      call cls1buf
      write(lun,'(a13,a4)') '%Feap-Macro  ',macro
c...  choose new color
      icolold=0
      return
      end
c
      subroutine gspst(string,l)
c-------------------------------------------------
c
C      Purpose: fill up buffer and clear if buffer is full
c
c      Input:
c       string - string
c       l      - position
c
c-------------------------------------------------
      USE hpgl1
      implicit double precision(a-h,o-z)
      character*(*) string
      character*132 gsbuff
      common /cpsbuf/ gsbuff
      common /ipsbuf/ ibuff,ibuf(481)
c.... buffer full (new line)
      if (l+ibuff.gt.132) then
         write(lun,'(a)') gsbuff(1:ibuff)
         ibuff=0
      end if
      gsbuff(ibuff+1:) = string(1:l)
      ibuff = ibuff+l
      return
      entry clsbuf
c.... write rest of buffer and final text
c     if (ibuff.gt.1) then
c       write(lun,'(a)') gsbuff(1:ibuff)
        write(lun,'(a)') 'stroke showpage '
c       ibuff=0
c     end if
      return
      entry cls1buf
c.... write rest of buffer to allow new text in psfile
      if (ibuff.gt.1) then
         write(lun,'(a)') gsbuff(1:ibuff)
         ibuff=0
      end if
      return
      end
c
      subroutine psfill(x,y,npts,nfill)
c-----------------------------------------------------------
c
c      Purpose: plot filled polygon
c
c      Input:
c      x(*)   - x-values of points
c      y(*)   - y-values of points
c      npts   - Number of points
c      nfill  - >0, fill polygon
c
c-----------------------------------------------------------
      implicit double precision(a-h,o-z)
      character*14 coord
      common /psplot1/ scl,xoff,yoff
      common /ipsbuf/  ibuff,ibuf(481)
      dimension x(*),y(*)
      data m1/-1/
      data ysl/480./
      ymx=y(1)
      ymn=y(1)
      do 10 i=2,npts
         ymx=max(ymx,y(i))
         ymn=min(ymn,y(i))
10    continue
      nscan = nint((ymx-ymn)*ysl)+2
      if(nfill.gt.0) call gspst('v ',2)
      do 20 i=1,nscan
         ibuf(i)=m1
20    continue
c.... move to
      xpp = scl*x(1)+xoff
      ypp = scl*y(1)+yoff
      write(coord,'(f5.1,1x,f5.1,1x,a2)')  xpp,ypp,'m '
      call gspst(coord,14)
c.... draw polygon
      do 30 i=2,npts
          xpp = scl*x(i)+xoff
          ypp = scl*y(i)+yoff
          write(coord,'(f5.1,1x,f5.1,1x,a2)')  xpp,ypp,'l '
          call gspst(coord,14)
30    continue
      if(nfill.gt.0) call gspst('f ',2)
      return
      end
C
      BLOCK DATA GSFONT
c------------------------------------------------------------
c
c     Purpose: this data is the font  for the hpgl text
c
c------------------------------------------------------------
      INTEGER*2 BXY
      COMMON/GCFONT/INDX(96),BXY(1556)
      DATA INDX(  1),INDX(  2),INDX(  3),INDX(  4)/   1,  17,  27,  47/
      DATA INDX(  5),INDX(  6),INDX(  7),INDX(  8)/  73, 100, 125, 138/
      DATA INDX(  9),INDX( 10),INDX( 11),INDX( 12)/ 147, 156, 171, 181/
      DATA INDX( 13),INDX( 14),INDX( 15),INDX( 16)/ 194, 199, 210, 215/
      DATA INDX( 17),INDX( 18),INDX( 19),INDX( 20)/ 236, 248, 269, 297/
      DATA INDX( 21),INDX( 22),INDX( 23),INDX( 24)/ 306, 325, 346, 359/
      DATA INDX( 25),INDX( 26),INDX( 27),INDX( 28)/ 395, 416, 438, 462/
      DATA INDX( 29),INDX( 30),INDX( 31),INDX( 32)/ 469, 479, 486, 514/
      DATA INDX( 33),INDX( 34),INDX( 35),INDX( 36)/ 543, 560, 587, 604/
      DATA INDX( 37),INDX( 38),INDX( 39),INDX( 40)/ 621, 635, 646, 667/
      DATA INDX( 41),INDX( 42),INDX( 43),INDX( 44)/ 681, 696, 712, 725/
      DATA INDX( 45),INDX( 46),INDX( 47),INDX( 48)/ 732, 747, 758, 777/
      DATA INDX( 49),INDX( 50),INDX( 51),INDX( 52)/ 791, 815, 834, 859/
      DATA INDX( 53),INDX( 54),INDX( 55),INDX( 56)/ 869, 882, 893, 909/
      DATA INDX( 57),INDX( 58),INDX( 59),INDX( 60)/ 926, 942, 955, 964/
      DATA INDX( 61),INDX( 62),INDX( 63),INDX( 64)/ 969, 978, 985, 990/
      DATA INDX( 65),INDX( 66),INDX( 67),INDX( 68)/1003,1029,1050,1067/
      DATA INDX( 69),INDX( 70),INDX( 71),INDX( 72)/1089,1108,1124,1152/
      DATA INDX( 73),INDX( 74),INDX( 75),INDX( 76)/1167,1190,1214,1228/
      DATA INDX( 77),INDX( 78),INDX( 79),INDX( 80)/1240,1261,1276,1295/
      DATA INDX( 81),INDX( 82),INDX( 83),INDX( 84)/1317,1339,1352,1375/
      DATA INDX( 85),INDX( 86),INDX( 87),INDX( 88)/1391,1407,1418,1442/
      DATA INDX( 89),INDX( 90),INDX( 91),INDX( 92)/1451,1473,1482,1501/
      DATA INDX( 93),INDX( 94),INDX( 95),INDX( 96)/1506,1525,1538,1555/
      DATA BXY(   1),BXY(   2),BXY(   3),BXY(   4)/-64,  3,  8,  3/
      DATA BXY(   5),BXY(   6),BXY(   7),BXY(   8)/  3,-64,  3,  0/
      DATA BXY(   9),BXY(  10),BXY(  11),BXY(  12)/  2,  0,  2,  1/
      DATA BXY(  13),BXY(  14),BXY(  15),BXY(  16)/  3,  1,  3,  0/
      DATA BXY(  17),BXY(  18),BXY(  19),BXY(  20)/-64,  2,  8,  2/
      DATA BXY(  21),BXY(  22),BXY(  23),BXY(  24)/  6,-64,  4,  8/
      DATA BXY(  25),BXY(  26),BXY(  27),BXY(  28)/  4,  6,-64,  2/
      DATA BXY(  29),BXY(  30),BXY(  31),BXY(  32)/  8,  2,  0,-64/
      DATA BXY(  33),BXY(  34),BXY(  35),BXY(  36)/  4,  8,  4,  0/
      DATA BXY(  37),BXY(  38),BXY(  39),BXY(  40)/-64,  6,  5,  0/
      DATA BXY(  41),BXY(  42),BXY(  43),BXY(  44)/  5,-64,  0,  3/
      DATA BXY(  45),BXY(  46),BXY(  47),BXY(  48)/  6,  3,-64,  6/
      DATA BXY(  49),BXY(  50),BXY(  51),BXY(  52)/  7,  1,  7,  0/
      DATA BXY(  53),BXY(  54),BXY(  55),BXY(  56)/  6,  0,  5,  1/
      DATA BXY(  57),BXY(  58),BXY(  59),BXY(  60)/  4,  5,  4,  6/
      DATA BXY(  61),BXY(  62),BXY(  63),BXY(  64)/  3,  6,  2,  5/
      DATA BXY(  65),BXY(  66),BXY(  67),BXY(  68)/  1,  0,  1,-64/
      DATA BXY(  69),BXY(  70),BXY(  71),BXY(  72)/  3,  8,  3,  0/
      DATA BXY(  73),BXY(  74),BXY(  75),BXY(  76)/-64,  1,  8,  0/
      DATA BXY(  77),BXY(  78),BXY(  79),BXY(  80)/  7,  1,  6,  2/
      DATA BXY(  81),BXY(  82),BXY(  83),BXY(  84)/  7,  1,  8,-64/
      DATA BXY(  85),BXY(  86),BXY(  87),BXY(  88)/  6,  7,  0,  1/
      DATA BXY(  89),BXY(  90),BXY(  91),BXY(  92)/-64,  5,  2,  4/
      DATA BXY(  93),BXY(  94),BXY(  95),BXY(  96)/  1,  5,  0,  6/
      DATA BXY(  97),BXY(  98),BXY(  99),BXY( 100)/  1,  5,  2,-64/
      DATA BXY( 101),BXY( 102),BXY( 103),BXY( 104)/  6,  3,  3,  0/
      DATA BXY( 105),BXY( 106),BXY( 107),BXY( 108)/  1,  0,  0,  1/
      DATA BXY( 109),BXY( 110),BXY( 111),BXY( 112)/  0,  2,  4,  6/
      DATA BXY( 113),BXY( 114),BXY( 115),BXY( 116)/  4,  7,  3,  8/
      DATA BXY( 117),BXY( 118),BXY( 119),BXY( 120)/  1,  8,  0,  7/
      DATA BXY( 121),BXY( 122),BXY( 123),BXY( 124)/  0,  6,  6,  0/
      DATA BXY( 125),BXY( 126),BXY( 127),BXY( 128)/-64,  4,  7,  3/
      DATA BXY( 129),BXY( 130),BXY( 131),BXY( 132)/  7,  3,  8,  4/
      DATA BXY( 133),BXY( 134),BXY( 135),BXY( 136)/  8,  4,  7,  2/
      DATA BXY( 137),BXY( 138),BXY( 139),BXY( 140)/  5,-64,  4,  8/
      DATA BXY( 141),BXY( 142),BXY( 143),BXY( 144)/  2,  6,  2,  2/
      DATA BXY( 145),BXY( 146),BXY( 147),BXY( 148)/  4,  0,-64,  2/
      DATA BXY( 149),BXY( 150),BXY( 151),BXY( 152)/  8,  4,  6,  4/
      DATA BXY( 153),BXY( 154),BXY( 155),BXY( 156)/  2,  2,  0,-64/
      DATA BXY( 157),BXY( 158),BXY( 159),BXY( 160)/  1,  2,  5,  6/
      DATA BXY( 161),BXY( 162),BXY( 163),BXY( 164)/-64,  3,  7,  3/
      DATA BXY( 165),BXY( 166),BXY( 167),BXY( 168)/  1,-64,  1,  6/
      DATA BXY( 169),BXY( 170),BXY( 171),BXY( 172)/  5,  2,-64,  3/
      DATA BXY( 173),BXY( 174),BXY( 175),BXY( 176)/  7,  3,  1,-64/
      DATA BXY( 177),BXY( 178),BXY( 179),BXY( 180)/  0,  4,  6,  4/
      DATA BXY( 181),BXY( 182),BXY( 183),BXY( 184)/-64,  3,  0,  2/
      DATA BXY( 185),BXY( 186),BXY( 187),BXY( 188)/  0,  2,  1,  3/
      DATA BXY( 189),BXY( 190),BXY( 191),BXY( 192)/  1,  3,  0,  1/
      DATA BXY( 193),BXY( 194),BXY( 195),BXY( 196)/ -2,-64,  0,  4/
      DATA BXY( 197),BXY( 198),BXY( 199),BXY( 200)/  6,  4,-64,  3/
      DATA BXY( 201),BXY( 202),BXY( 203),BXY( 204)/  0,  2,  0,  2/
      DATA BXY( 205),BXY( 206),BXY( 207),BXY( 208)/  1,  3,  1,  3/
      DATA BXY( 209),BXY( 210),BXY( 211),BXY( 212)/  0,-64,  0,  1/
      DATA BXY( 213),BXY( 214),BXY( 215),BXY( 216)/  6,  7,-64,  6/
      DATA BXY( 217),BXY( 218),BXY( 219),BXY( 220)/  7,  6,  1,  5/
      DATA BXY( 221),BXY( 222),BXY( 223),BXY( 224)/  0,  1,  0,  0/
      DATA BXY( 225),BXY( 226),BXY( 227),BXY( 228)/  1,  0,  7,  1/
      DATA BXY( 229),BXY( 230),BXY( 231),BXY( 232)/  8,  5,  8,  6/
      DATA BXY( 233),BXY( 234),BXY( 235),BXY( 236)/  7,  0,  1,-64/
      DATA BXY( 237),BXY( 238),BXY( 239),BXY( 240)/  1,  6,  3,  8/
      DATA BXY( 241),BXY( 242),BXY( 243),BXY( 244)/  3,  0,-64,  1/
      DATA BXY( 245),BXY( 246),BXY( 247),BXY( 248)/  0,  5,  0,-64/
      DATA BXY( 249),BXY( 250),BXY( 251),BXY( 252)/  0,  7,  1,  8/
      DATA BXY( 253),BXY( 254),BXY( 255),BXY( 256)/  5,  8,  6,  7/
      DATA BXY( 257),BXY( 258),BXY( 259),BXY( 260)/  6,  6,  4,  4/
      DATA BXY( 261),BXY( 262),BXY( 263),BXY( 264)/  2,  4,  0,  2/
      DATA BXY( 265),BXY( 266),BXY( 267),BXY( 268)/  0,  0,  6,  0/
      DATA BXY( 269),BXY( 270),BXY( 271),BXY( 272)/-64,  0,  7,  1/
      DATA BXY( 273),BXY( 274),BXY( 275),BXY( 276)/  8,  5,  8,  6/
      DATA BXY( 277),BXY( 278),BXY( 279),BXY( 280)/  7,  6,  5,  5/
      DATA BXY( 281),BXY( 282),BXY( 283),BXY( 284)/  4,  1,  4,-64/
      DATA BXY( 285),BXY( 286),BXY( 287),BXY( 288)/  5,  4,  6,  3/
      DATA BXY( 289),BXY( 290),BXY( 291),BXY( 292)/  6,  1,  5,  0/
      DATA BXY( 293),BXY( 294),BXY( 295),BXY( 296)/  1,  0,  0,  1/
      DATA BXY( 297),BXY( 298),BXY( 299),BXY( 300)/-64,  5,  0,  5/
      DATA BXY( 301),BXY( 302),BXY( 303),BXY( 304)/  8,  0,  3,  6/
      DATA BXY( 305),BXY( 306),BXY( 307),BXY( 308)/  3,-64,  0,  1/
      DATA BXY( 309),BXY( 310),BXY( 311),BXY( 312)/  1,  0,  4,  0/
      DATA BXY( 313),BXY( 314),BXY( 315),BXY( 316)/  6,  2,  6,  3/
      DATA BXY( 317),BXY( 318),BXY( 319),BXY( 320)/  4,  5,  0,  5/
      DATA BXY( 321),BXY( 322),BXY( 323),BXY( 324)/  0,  8,  6,  8/
      DATA BXY( 325),BXY( 326),BXY( 327),BXY( 328)/-64,  5,  8,  2/
      DATA BXY( 329),BXY( 330),BXY( 331),BXY( 332)/  8,  0,  6,  0/
      DATA BXY( 333),BXY( 334),BXY( 335),BXY( 336)/  1,  1,  0,  5/
      DATA BXY( 337),BXY( 338),BXY( 339),BXY( 340)/  0,  6,  1,  6/
      DATA BXY( 341),BXY( 342),BXY( 343),BXY( 344)/  3,  5,  4,  0/
      DATA BXY( 345),BXY( 346),BXY( 347),BXY( 348)/  4,-64,  0,  7/
      DATA BXY( 349),BXY( 350),BXY( 351),BXY( 352)/  0,  8,  6,  8/
      DATA BXY( 353),BXY( 354),BXY( 355),BXY( 356)/  6,  7,  2,  3/
      DATA BXY( 357),BXY( 358),BXY( 359),BXY( 360)/  2,  0,-64,  6/
      DATA BXY( 361),BXY( 362),BXY( 363),BXY( 364)/  7,  5,  8,  1/
      DATA BXY( 365),BXY( 366),BXY( 367),BXY( 368)/  8,  0,  7,  0/
      DATA BXY( 369),BXY( 370),BXY( 371),BXY( 372)/  5,  1,  4,  5/
      DATA BXY( 373),BXY( 374),BXY( 375),BXY( 376)/  4,  6,  3,  6/
      DATA BXY( 377),BXY( 378),BXY( 379),BXY( 380)/  1,  5,  0,  1/
      DATA BXY( 381),BXY( 382),BXY( 383),BXY( 384)/  0,  0,  1,  0/
      DATA BXY( 385),BXY( 386),BXY( 387),BXY( 388)/  3,  1,  4,-64/
      DATA BXY( 389),BXY( 390),BXY( 391),BXY( 392)/  5,  4,  6,  5/
      DATA BXY( 393),BXY( 394),BXY( 395),BXY( 396)/  6,  7,-64,  1/
      DATA BXY( 397),BXY( 398),BXY( 399),BXY( 400)/  0,  4,  0,  6/
      DATA BXY( 401),BXY( 402),BXY( 403),BXY( 404)/  2,  6,  7,  5/
      DATA BXY( 405),BXY( 406),BXY( 407),BXY( 408)/  8,  1,  8,  0/
      DATA BXY( 409),BXY( 410),BXY( 411),BXY( 412)/  7,  0,  5,  1/
      DATA BXY( 413),BXY( 414),BXY( 415),BXY( 416)/  4,  6,  4,-64/
      DATA BXY( 417),BXY( 418),BXY( 419),BXY( 420)/  3,  4,  2,  4/
      DATA BXY( 421),BXY( 422),BXY( 423),BXY( 424)/  2,  5,  3,  5/
      DATA BXY( 425),BXY( 426),BXY( 427),BXY( 428)/  3,  4,-64,  3/
      DATA BXY( 429),BXY( 430),BXY( 431),BXY( 432)/  0,  2,  0,  2/
      DATA BXY( 433),BXY( 434),BXY( 435),BXY( 436)/  1,  3,  1,  3/
      DATA BXY( 437),BXY( 438),BXY( 439),BXY( 440)/  0,-64,  3,  4/
      DATA BXY( 441),BXY( 442),BXY( 443),BXY( 444)/  2,  4,  2,  5/
      DATA BXY( 445),BXY( 446),BXY( 447),BXY( 448)/  3,  5,  3,  4/
      DATA BXY( 449),BXY( 450),BXY( 451),BXY( 452)/-64,  3,  0,  2/
      DATA BXY( 453),BXY( 454),BXY( 455),BXY( 456)/  0,  2,  1,  3/
      DATA BXY( 457),BXY( 458),BXY( 459),BXY( 460)/  1,  3,  0,  1/
      DATA BXY( 461),BXY( 462),BXY( 463),BXY( 464)/ -2,-64,  4,  8/
      DATA BXY( 465),BXY( 466),BXY( 467),BXY( 468)/  0,  4,  4,  0/
      DATA BXY( 469),BXY( 470),BXY( 471),BXY( 472)/-64,  5,  5,  1/
      DATA BXY( 473),BXY( 474),BXY( 475),BXY( 476)/  5,-64,  1,  3/
      DATA BXY( 477),BXY( 478),BXY( 479),BXY( 480)/  5,  3,-64,  2/
      DATA BXY( 481),BXY( 482),BXY( 483),BXY( 484)/  8,  6,  4,  2/
      DATA BXY( 485),BXY( 486),BXY( 487),BXY( 488)/  0,-64,  0,  7/
      DATA BXY( 489),BXY( 490),BXY( 491),BXY( 492)/  1,  8,  5,  8/
      DATA BXY( 493),BXY( 494),BXY( 495),BXY( 496)/  6,  7,  6,  5/
      DATA BXY( 497),BXY( 498),BXY( 499),BXY( 500)/  5,  4,  3,  4/
      DATA BXY( 501),BXY( 502),BXY( 503),BXY( 504)/  3,  3,-64,  3/
      DATA BXY( 505),BXY( 506),BXY( 507),BXY( 508)/  0,  2,  0,  2/
      DATA BXY( 509),BXY( 510),BXY( 511),BXY( 512)/  1,  3,  1,  3/
      DATA BXY( 513),BXY( 514),BXY( 515),BXY( 516)/  0,-64,  4,  3/
      DATA BXY( 517),BXY( 518),BXY( 519),BXY( 520)/  4,  5,  3,  5/
      DATA BXY( 521),BXY( 522),BXY( 523),BXY( 524)/  2,  4,  2,  3/
      DATA BXY( 525),BXY( 526),BXY( 527),BXY( 528)/  5,  3,  6,  4/
      DATA BXY( 529),BXY( 530),BXY( 531),BXY( 532)/  6,  7,  5,  8/
      DATA BXY( 533),BXY( 534),BXY( 535),BXY( 536)/  2,  8,  0,  6/
      DATA BXY( 537),BXY( 538),BXY( 539),BXY( 540)/  0,  2,  2,  0/
      DATA BXY( 541),BXY( 542),BXY( 543),BXY( 544)/  5,  0,  0,  0/
      DATA BXY( 545),BXY( 546),BXY( 547),BXY( 548)/  0,  6,  2,  8/
      DATA BXY( 549),BXY( 550),BXY( 551),BXY( 552)/  4,  8,  6,  6/
      DATA BXY( 553),BXY( 554),BXY( 555),BXY( 556)/  6,  0,-64,  0/
      DATA BXY( 557),BXY( 558),BXY( 559),BXY( 560)/  3,  6,  3,  0/
      DATA BXY( 561),BXY( 562),BXY( 563),BXY( 564)/  0,  5,  0,  6/
      DATA BXY( 565),BXY( 566),BXY( 567),BXY( 568)/  1,  6,  3,  5/
      DATA BXY( 569),BXY( 570),BXY( 571),BXY( 572)/  4,  1,  4,-64/
      DATA BXY( 573),BXY( 574),BXY( 575),BXY( 576)/  5,  4,  6,  5/
      DATA BXY( 577),BXY( 578),BXY( 579),BXY( 580)/  6,  7,  5,  8/
      DATA BXY( 581),BXY( 582),BXY( 583),BXY( 584)/  0,  8,  1,  8/
      DATA BXY( 585),BXY( 586),BXY( 587),BXY( 588)/  1,  0,-64,  6/
      DATA BXY( 589),BXY( 590),BXY( 591),BXY( 592)/  7,  5,  8,  2/
      DATA BXY( 593),BXY( 594),BXY( 595),BXY( 596)/  8,  0,  6,  0/
      DATA BXY( 597),BXY( 598),BXY( 599),BXY( 600)/  2,  2,  0,  5/
      DATA BXY( 601),BXY( 602),BXY( 603),BXY( 604)/  0,  6,  1,  0/
      DATA BXY( 605),BXY( 606),BXY( 607),BXY( 608)/  0,  4,  0,  6/
      DATA BXY( 609),BXY( 610),BXY( 611),BXY( 612)/  2,  6,  6,  4/
      DATA BXY( 613),BXY( 614),BXY( 615),BXY( 616)/  8,  0,  8,-64/
      DATA BXY( 617),BXY( 618),BXY( 619),BXY( 620)/  1,  8,  1,  0/
      DATA BXY( 621),BXY( 622),BXY( 623),BXY( 624)/-64,  6,  0,  0/
      DATA BXY( 625),BXY( 626),BXY( 627),BXY( 628)/  0,  0,  8,  6/
      DATA BXY( 629),BXY( 630),BXY( 631),BXY( 632)/  8,-64,  3,  4/
      DATA BXY( 633),BXY( 634),BXY( 635),BXY( 636)/  0,  4,  0,  0/
      DATA BXY( 637),BXY( 638),BXY( 639),BXY( 640)/  0,  8,  6,  8/
      DATA BXY( 641),BXY( 642),BXY( 643),BXY( 644)/-64,  3,  4,  0/
      DATA BXY( 645),BXY( 646),BXY( 647),BXY( 648)/  4,-64,  6,  7/
      DATA BXY( 649),BXY( 650),BXY( 651),BXY( 652)/  5,  8,  2,  8/
      DATA BXY( 653),BXY( 654),BXY( 655),BXY( 656)/  0,  6,  0,  2/
      DATA BXY( 657),BXY( 658),BXY( 659),BXY( 660)/  2,  0,  5,  0/
      DATA BXY( 661),BXY( 662),BXY( 663),BXY( 664)/  6,  1,  6,  3/
      DATA BXY( 665),BXY( 666),BXY( 667),BXY( 668)/  3,  3,  0,  0/
      DATA BXY( 669),BXY( 670),BXY( 671),BXY( 672)/  0,  8,-64,  0/
      DATA BXY( 673),BXY( 674),BXY( 675),BXY( 676)/  4,  6,  4,-64/
      DATA BXY( 677),BXY( 678),BXY( 679),BXY( 680)/  6,  8,  6,  0/
      DATA BXY( 681),BXY( 682),BXY( 683),BXY( 684)/-64,  1,  0,  5/
      DATA BXY( 685),BXY( 686),BXY( 687),BXY( 688)/  0,-64,  3,  0/
      DATA BXY( 689),BXY( 690),BXY( 691),BXY( 692)/  3,  8,-64,  1/
      DATA BXY( 693),BXY( 694),BXY( 695),BXY( 696)/  8,  5,  8,-64/
      DATA BXY( 697),BXY( 698),BXY( 699),BXY( 700)/  0,  1,  1,  0/
      DATA BXY( 701),BXY( 702),BXY( 703),BXY( 704)/  3,  0,  4,  1/
      DATA BXY( 705),BXY( 706),BXY( 707),BXY( 708)/  4,  8,-64,  2/
      DATA BXY( 709),BXY( 710),BXY( 711),BXY( 712)/  8,  6,  8,  0/
      DATA BXY( 713),BXY( 714),BXY( 715),BXY( 716)/  0,  0,  8,-64/
      DATA BXY( 717),BXY( 718),BXY( 719),BXY( 720)/  6,  8,  0,  2/
      DATA BXY( 721),BXY( 722),BXY( 723),BXY( 724)/  2,  4,  6,  0/
      DATA BXY( 725),BXY( 726),BXY( 727),BXY( 728)/-64,  0,  8,  0/
      DATA BXY( 729),BXY( 730),BXY( 731),BXY( 732)/  0,  6,  0,  0/
      DATA BXY( 733),BXY( 734),BXY( 735),BXY( 736)/  0,  0,  8,  3/
      DATA BXY( 737),BXY( 738),BXY( 739),BXY( 740)/  5,  3,  4,-64/
      DATA BXY( 741),BXY( 742),BXY( 743),BXY( 744)/  3,  5,  6,  8/
      DATA BXY( 745),BXY( 746),BXY( 747),BXY( 748)/  6,  0,  0,  0/
      DATA BXY( 749),BXY( 750),BXY( 751),BXY( 752)/  0,  8,  6,  2/
      DATA BXY( 753),BXY( 754),BXY( 755),BXY( 756)/-64,  6,  8,  6/
      DATA BXY( 757),BXY( 758),BXY( 759),BXY( 760)/  0,-64,  6,  2/
      DATA BXY( 761),BXY( 762),BXY( 763),BXY( 764)/  6,  6,  4,  8/
      DATA BXY( 765),BXY( 766),BXY( 767),BXY( 768)/  2,  8,  0,  6/
      DATA BXY( 769),BXY( 770),BXY( 771),BXY( 772)/  0,  2,  2,  0/
      DATA BXY( 773),BXY( 774),BXY( 775),BXY( 776)/  4,  0,  6,  2/
      DATA BXY( 777),BXY( 778),BXY( 779),BXY( 780)/  0,  0,  0,  8/
      DATA BXY( 781),BXY( 782),BXY( 783),BXY( 784)/  5,  8,  6,  7/
      DATA BXY( 785),BXY( 786),BXY( 787),BXY( 788)/  6,  5,  5,  4/
      DATA BXY( 789),BXY( 790),BXY( 791),BXY( 792)/  0,  4,-64,  6/
      DATA BXY( 793),BXY( 794),BXY( 795),BXY( 796)/  2,  6,  6,  4/
      DATA BXY( 797),BXY( 798),BXY( 799),BXY( 800)/  8,  2,  8,  0/
      DATA BXY( 801),BXY( 802),BXY( 803),BXY( 804)/  6,  0,  2,  2/
      DATA BXY( 805),BXY( 806),BXY( 807),BXY( 808)/  0,  4,  0,  6/
      DATA BXY( 809),BXY( 810),BXY( 811),BXY( 812)/  2,-64,  3,  3/
      DATA BXY( 813),BXY( 814),BXY( 815),BXY( 816)/  6,  0,  0,  0/
      DATA BXY( 817),BXY( 818),BXY( 819),BXY( 820)/  0,  8,  5,  8/
      DATA BXY( 821),BXY( 822),BXY( 823),BXY( 824)/  6,  7,  6,  5/
      DATA BXY( 825),BXY( 826),BXY( 827),BXY( 828)/  5,  4,  0,  4/
      DATA BXY( 829),BXY( 830),BXY( 831),BXY( 832)/-64,  2,  4,  6/
      DATA BXY( 833),BXY( 834),BXY( 835),BXY( 836)/  0,-64,  6,  7/
      DATA BXY( 837),BXY( 838),BXY( 839),BXY( 840)/  5,  8,  1,  8/
      DATA BXY( 841),BXY( 842),BXY( 843),BXY( 844)/  0,  7,  0,  5/
      DATA BXY( 845),BXY( 846),BXY( 847),BXY( 848)/  1,  4,  5,  4/
      DATA BXY( 849),BXY( 850),BXY( 851),BXY( 852)/  6,  3,  6,  1/
      DATA BXY( 853),BXY( 854),BXY( 855),BXY( 856)/  5,  0,  1,  0/
      DATA BXY( 857),BXY( 858),BXY( 859),BXY( 860)/  0,  1,-64,  0/
      DATA BXY( 861),BXY( 862),BXY( 863),BXY( 864)/  8,  6,  8,-64/
      DATA BXY( 865),BXY( 866),BXY( 867),BXY( 868)/  3,  0,  3,  8/
      DATA BXY( 869),BXY( 870),BXY( 871),BXY( 872)/-64,  6,  8,  6/
      DATA BXY( 873),BXY( 874),BXY( 875),BXY( 876)/  1,  5,  0,  1/
      DATA BXY( 877),BXY( 878),BXY( 879),BXY( 880)/  0,  0,  1,  0/
      DATA BXY( 881),BXY( 882),BXY( 883),BXY( 884)/  8,-64,  0,  8/
      DATA BXY( 885),BXY( 886),BXY( 887),BXY( 888)/  0,  6,  3,  0/
      DATA BXY( 889),BXY( 890),BXY( 891),BXY( 892)/  6,  6,  6,  8/
      DATA BXY( 893),BXY( 894),BXY( 895),BXY( 896)/-64,  0,  8,  0/
      DATA BXY( 897),BXY( 898),BXY( 899),BXY( 900)/  0,  3,  3,  3/
      DATA BXY( 901),BXY( 902),BXY( 903),BXY( 904)/  4,-64,  3,  3/
      DATA BXY( 905),BXY( 906),BXY( 907),BXY( 908)/  6,  0,  6,  8/
      DATA BXY( 909),BXY( 910),BXY( 911),BXY( 912)/  0,  0,  0,  1/
      DATA BXY( 913),BXY( 914),BXY( 915),BXY( 916)/  6,  7,  6,  8/
      DATA BXY( 917),BXY( 918),BXY( 919),BXY( 920)/-64,  0,  8,  0/
      DATA BXY( 921),BXY( 922),BXY( 923),BXY( 924)/  7,  6,  1,  6/
      DATA BXY( 925),BXY( 926),BXY( 927),BXY( 928)/  0,-64,  0,  8/
      DATA BXY( 929),BXY( 930),BXY( 931),BXY( 932)/  0,  7,  3,  4/
      DATA BXY( 933),BXY( 934),BXY( 935),BXY( 936)/  6,  7,  6,  8/
      DATA BXY( 937),BXY( 938),BXY( 939),BXY( 940)/-64,  3,  4,  3/
      DATA BXY( 941),BXY( 942),BXY( 943),BXY( 944)/  0,-64,  0,  8/
      DATA BXY( 945),BXY( 946),BXY( 947),BXY( 948)/  6,  8,  6,  7/
      DATA BXY( 949),BXY( 950),BXY( 951),BXY( 952)/  0,  1,  0,  0/
      DATA BXY( 953),BXY( 954),BXY( 955),BXY( 956)/  6,  0,-64,  4/
      DATA BXY( 957),BXY( 958),BXY( 959),BXY( 960)/  8,  2,  8,  2/
      DATA BXY( 961),BXY( 962),BXY( 963),BXY( 964)/  0,  4,  0,-64/
      DATA BXY( 965),BXY( 966),BXY( 967),BXY( 968)/  0,  7,  6,  1/
      DATA BXY( 969),BXY( 970),BXY( 971),BXY( 972)/-64,  3,  8,  5/
      DATA BXY( 973),BXY( 974),BXY( 975),BXY( 976)/  8,  5,  0,  3/
      DATA BXY( 977),BXY( 978),BXY( 979),BXY( 980)/  0,-64,  0,  5/
      DATA BXY( 981),BXY( 982),BXY( 983),BXY( 984)/  3,  8,  6,  5/
      DATA BXY( 985),BXY( 986),BXY( 987),BXY( 988)/-64,  0, -2,  6/
      DATA BXY( 989),BXY( 990),BXY( 991),BXY( 992)/ -2,-64,  2,  7/
      DATA BXY( 993),BXY( 994),BXY( 995),BXY( 996)/  3,  7,  3,  8/
      DATA BXY( 997),BXY( 998),BXY( 999),BXY(1000)/  2,  8,  2,  7/
      DATA BXY(1001),BXY(1002),BXY(1003),BXY(1004)/  4,  5,-64,  5/
      DATA BXY(1005),BXY(1006),BXY(1007),BXY(1008)/  3,  1,  3,  0/
      DATA BXY(1009),BXY(1010),BXY(1011),BXY(1012)/  2,  0,  1,  1/
      DATA BXY(1013),BXY(1014),BXY(1015),BXY(1016)/  0,  4,  0,  5/
      DATA BXY(1017),BXY(1018),BXY(1019),BXY(1020)/  1,-64,  1,  5/
      DATA BXY(1021),BXY(1022),BXY(1023),BXY(1024)/  4,  5,  5,  4/
      DATA BXY(1025),BXY(1026),BXY(1027),BXY(1028)/  5,  1,  6,  0/
      DATA BXY(1029),BXY(1030),BXY(1031),BXY(1032)/  0,  0,  0,  8/
      DATA BXY(1033),BXY(1034),BXY(1035),BXY(1036)/-64,  0,  3,  2/
      DATA BXY(1037),BXY(1038),BXY(1039),BXY(1040)/  5,  5,  5,  6/
      DATA BXY(1041),BXY(1042),BXY(1043),BXY(1044)/  4,  6,  1,  5/
      DATA BXY(1045),BXY(1046),BXY(1047),BXY(1048)/  0,  2,  0,  0/
      DATA BXY(1049),BXY(1050),BXY(1051),BXY(1052)/  2,-64,  5,  4/
      DATA BXY(1053),BXY(1054),BXY(1055),BXY(1056)/  4,  5,  1,  5/
      DATA BXY(1057),BXY(1058),BXY(1059),BXY(1060)/  0,  4,  0,  1/
      DATA BXY(1061),BXY(1062),BXY(1063),BXY(1064)/  1,  0,  4,  0/
      DATA BXY(1065),BXY(1066),BXY(1067),BXY(1068)/  5,  1,-64,  5/
      DATA BXY(1069),BXY(1070),BXY(1071),BXY(1072)/  8,  5,  0,-64/
      DATA BXY(1073),BXY(1074),BXY(1075),BXY(1076)/  5,  2,  3,  0/
      DATA BXY(1077),BXY(1078),BXY(1079),BXY(1080)/  1,  0,  0,  1/
      DATA BXY(1081),BXY(1082),BXY(1083),BXY(1084)/  0,  4,  1,  5/
      DATA BXY(1085),BXY(1086),BXY(1087),BXY(1088)/  3,  5,  5,  3/
      DATA BXY(1089),BXY(1090),BXY(1091),BXY(1092)/-64,  0,  3,  6/
      DATA BXY(1093),BXY(1094),BXY(1095),BXY(1096)/  3,  6,  4,  5/
      DATA BXY(1097),BXY(1098),BXY(1099),BXY(1100)/  5,  1,  5,  0/
      DATA BXY(1101),BXY(1102),BXY(1103),BXY(1104)/  4,  0,  1,  1/
      DATA BXY(1105),BXY(1106),BXY(1107),BXY(1108)/  0,  5,  0,-64/
      DATA BXY(1109),BXY(1110),BXY(1111),BXY(1112)/  2,  0,  2,  7/
      DATA BXY(1113),BXY(1114),BXY(1115),BXY(1116)/  3,  8,  4,  8/
      DATA BXY(1117),BXY(1118),BXY(1119),BXY(1120)/  5,  7,-64,  0/
      DATA BXY(1121),BXY(1122),BXY(1123),BXY(1124)/  4,  4,  4,-64/
      DATA BXY(1125),BXY(1126),BXY(1127),BXY(1128)/  5,  2,  3,  0/
      DATA BXY(1129),BXY(1130),BXY(1131),BXY(1132)/  1,  0,  0,  1/
      DATA BXY(1133),BXY(1134),BXY(1135),BXY(1136)/  0,  4,  1,  5/
      DATA BXY(1137),BXY(1138),BXY(1139),BXY(1140)/  3,  5,  5,  3/
      DATA BXY(1141),BXY(1142),BXY(1143),BXY(1144)/-64,  5,  5,  5/
      DATA BXY(1145),BXY(1146),BXY(1147),BXY(1148)/ -3,  4, -4,  1/
      DATA BXY(1149),BXY(1150),BXY(1151),BXY(1152)/ -4,  0, -3,  0/
      DATA BXY(1153),BXY(1154),BXY(1155),BXY(1156)/  0,  0,  8,-64/
      DATA BXY(1157),BXY(1158),BXY(1159),BXY(1160)/  0,  3,  2,  5/
      DATA BXY(1161),BXY(1162),BXY(1163),BXY(1164)/  5,  5,  6,  4/
      DATA BXY(1165),BXY(1166),BXY(1167),BXY(1168)/  6,  0,-64,  2/
      DATA BXY(1169),BXY(1170),BXY(1171),BXY(1172)/  0,  4,  0,-64/
      DATA BXY(1173),BXY(1174),BXY(1175),BXY(1176)/  3,  0,  3,  5/
      DATA BXY(1177),BXY(1178),BXY(1179),BXY(1180)/  2,  5,-64,  2/
      DATA BXY(1181),BXY(1182),BXY(1183),BXY(1184)/  7,  3,  7,  3/
      DATA BXY(1185),BXY(1186),BXY(1187),BXY(1188)/  8,  2,  8,  2/
      DATA BXY(1189),BXY(1190),BXY(1191),BXY(1192)/  7,-64,  5,  7/
      DATA BXY(1193),BXY(1194),BXY(1195),BXY(1196)/  5,  8,  4,  8/
      DATA BXY(1197),BXY(1198),BXY(1199),BXY(1200)/  4,  7,  5,  7/
      DATA BXY(1201),BXY(1202),BXY(1203),BXY(1204)/-64,  4,  5,  5/
      DATA BXY(1205),BXY(1206),BXY(1207),BXY(1208)/  5,  5, -3,  4/
      DATA BXY(1209),BXY(1210),BXY(1211),BXY(1212)/ -4,  2, -4,  1/
      DATA BXY(1213),BXY(1214),BXY(1215),BXY(1216)/ -3,  0,  0,  0/
      DATA BXY(1217),BXY(1218),BXY(1219),BXY(1220)/  8,-64,  4,  5/
      DATA BXY(1221),BXY(1222),BXY(1223),BXY(1224)/  0,  1,-64,  2/
      DATA BXY(1225),BXY(1226),BXY(1227),BXY(1228)/  3,  5,  0,-64/
      DATA BXY(1229),BXY(1230),BXY(1231),BXY(1232)/  2,  0,  4,  0/
      DATA BXY(1233),BXY(1234),BXY(1235),BXY(1236)/-64,  3,  0,  3/
      DATA BXY(1237),BXY(1238),BXY(1239),BXY(1240)/  8,  2,  8,  0/
      DATA BXY(1241),BXY(1242),BXY(1243),BXY(1244)/  0,  0,  5,  2/
      DATA BXY(1245),BXY(1246),BXY(1247),BXY(1248)/  5,  3,  4,  3/
      DATA BXY(1249),BXY(1250),BXY(1251),BXY(1252)/  0,-64,  3,  4/
      DATA BXY(1253),BXY(1254),BXY(1255),BXY(1256)/  4,  5,  5,  5/
      DATA BXY(1257),BXY(1258),BXY(1259),BXY(1260)/  6,  4,  6,  0/
      DATA BXY(1261),BXY(1262),BXY(1263),BXY(1264)/  0,  0,  0,  5/
      DATA BXY(1265),BXY(1266),BXY(1267),BXY(1268)/-64,  0,  3,  2/
      DATA BXY(1269),BXY(1270),BXY(1271),BXY(1272)/  5,  4,  5,  5/
      DATA BXY(1273),BXY(1274),BXY(1275),BXY(1276)/  4,  5,  0,-64/
      DATA BXY(1277),BXY(1278),BXY(1279),BXY(1280)/  5,  4,  4,  5/
      DATA BXY(1281),BXY(1282),BXY(1283),BXY(1284)/  1,  5,  0,  4/
      DATA BXY(1285),BXY(1286),BXY(1287),BXY(1288)/  0,  1,  1,  0/
      DATA BXY(1289),BXY(1290),BXY(1291),BXY(1292)/  4,  0,  5,  1/
      DATA BXY(1293),BXY(1294),BXY(1295),BXY(1296)/  5,  4,-64,  0/
      DATA BXY(1297),BXY(1298),BXY(1299),BXY(1300)/ -4,  0,  5,-64/
      DATA BXY(1301),BXY(1302),BXY(1303),BXY(1304)/  0,  3,  2,  5/
      DATA BXY(1305),BXY(1306),BXY(1307),BXY(1308)/  4,  5,  5,  4/
      DATA BXY(1309),BXY(1310),BXY(1311),BXY(1312)/  5,  1,  4,  0/
      DATA BXY(1313),BXY(1314),BXY(1315),BXY(1316)/  2,  0,  0,  2/
      DATA BXY(1317),BXY(1318),BXY(1319),BXY(1320)/-64,  5, -4,  5/
      DATA BXY(1321),BXY(1322),BXY(1323),BXY(1324)/  5,-64,  5,  2/
      DATA BXY(1325),BXY(1326),BXY(1327),BXY(1328)/  3,  0,  1,  0/
      DATA BXY(1329),BXY(1330),BXY(1331),BXY(1332)/  0,  1,  0,  4/
      DATA BXY(1333),BXY(1334),BXY(1335),BXY(1336)/  1,  5,  3,  5/
      DATA BXY(1337),BXY(1338),BXY(1339),BXY(1340)/  5,  3,  0,  0/
      DATA BXY(1341),BXY(1342),BXY(1343),BXY(1344)/  0,  5,-64,  0/
      DATA BXY(1345),BXY(1346),BXY(1347),BXY(1348)/  3,  2,  5,  4/
      DATA BXY(1349),BXY(1350),BXY(1351),BXY(1352)/  5,  5,  4,-64/
      DATA BXY(1353),BXY(1354),BXY(1355),BXY(1356)/-64,-64,  0,  1/
      DATA BXY(1357),BXY(1358),BXY(1359),BXY(1360)/  1,  0,  4,  0/
      DATA BXY(1361),BXY(1362),BXY(1363),BXY(1364)/  5,  1,  4,  2/
      DATA BXY(1365),BXY(1366),BXY(1367),BXY(1368)/  1,  3,  0,  4/
      DATA BXY(1369),BXY(1370),BXY(1371),BXY(1372)/  1,  5,  4,  5/
      DATA BXY(1373),BXY(1374),BXY(1375),BXY(1376)/  5,  4,-64,  2/
      DATA BXY(1377),BXY(1378),BXY(1379),BXY(1380)/  8,  2,  1,  3/
      DATA BXY(1381),BXY(1382),BXY(1383),BXY(1384)/  0,  4,  0,  5/
      DATA BXY(1385),BXY(1386),BXY(1387),BXY(1388)/  1,-64,  0,  5/
      DATA BXY(1389),BXY(1390),BXY(1391),BXY(1392)/  4,  5,-64,  0/
      DATA BXY(1393),BXY(1394),BXY(1395),BXY(1396)/  5,  0,  1,  1/
      DATA BXY(1397),BXY(1398),BXY(1399),BXY(1400)/  0,  4,  0,  5/
      DATA BXY(1401),BXY(1402),BXY(1403),BXY(1404)/  1,-64,  5,  5/
      DATA BXY(1405),BXY(1406),BXY(1407),BXY(1408)/  5,  0,-64,  0/
      DATA BXY(1409),BXY(1410),BXY(1411),BXY(1412)/  5,  0,  3,  3/
      DATA BXY(1413),BXY(1414),BXY(1415),BXY(1416)/  0,  6,  3,  6/
      DATA BXY(1417),BXY(1418),BXY(1419),BXY(1420)/  5,-64,  0,  5/
      DATA BXY(1421),BXY(1422),BXY(1423),BXY(1424)/  0,  1,  1,  0/
      DATA BXY(1425),BXY(1426),BXY(1427),BXY(1428)/  2,  0,  3,  1/
      DATA BXY(1429),BXY(1430),BXY(1431),BXY(1432)/  3,  4,-64,  3/
      DATA BXY(1433),BXY(1434),BXY(1435),BXY(1436)/  1,  4,  0,  5/
      DATA BXY(1437),BXY(1438),BXY(1439),BXY(1440)/  0,  6,  1,  6/
      DATA BXY(1441),BXY(1442),BXY(1443),BXY(1444)/  5,  0,  0,  5/
      DATA BXY(1445),BXY(1446),BXY(1447),BXY(1448)/  5,-64,  0,  5/
      DATA BXY(1449),BXY(1450),BXY(1451),BXY(1452)/  5,  0,-64,  0/
      DATA BXY(1453),BXY(1454),BXY(1455),BXY(1456)/  5,  0,  1,  1/
      DATA BXY(1457),BXY(1458),BXY(1459),BXY(1460)/  0,  4,  0,  5/
      DATA BXY(1461),BXY(1462),BXY(1463),BXY(1464)/  1,-64,  5,  5/
      DATA BXY(1465),BXY(1466),BXY(1467),BXY(1468)/  5, -3,  4, -4/
      DATA BXY(1469),BXY(1470),BXY(1471),BXY(1472)/  1, -4,  0, -3/
      DATA BXY(1473),BXY(1474),BXY(1475),BXY(1476)/-64,  0,  5,  5/
      DATA BXY(1477),BXY(1478),BXY(1479),BXY(1480)/  5,  0,  0,  5/
      DATA BXY(1481),BXY(1482),BXY(1483),BXY(1484)/  0,-64,  3,  8/
      DATA BXY(1485),BXY(1486),BXY(1487),BXY(1488)/  2,  8,  1,  7/
      DATA BXY(1489),BXY(1490),BXY(1491),BXY(1492)/  1,  5,  0,  4/
      DATA BXY(1493),BXY(1494),BXY(1495),BXY(1496)/  1,  3,  1,  1/
      DATA BXY(1497),BXY(1498),BXY(1499),BXY(1500)/  2,  0,  3,  0/
      DATA BXY(1501),BXY(1502),BXY(1503),BXY(1504)/-64,  3,  8,  3/
      DATA BXY(1505),BXY(1506),BXY(1507),BXY(1508)/  0,-64,  3,  8/
      DATA BXY(1509),BXY(1510),BXY(1511),BXY(1512)/  4,  8,  5,  7/
      DATA BXY(1513),BXY(1514),BXY(1515),BXY(1516)/  5,  5,  6,  4/
      DATA BXY(1517),BXY(1518),BXY(1519),BXY(1520)/  5,  3,  5,  1/
      DATA BXY(1521),BXY(1522),BXY(1523),BXY(1524)/  4,  0,  3,  0/
      DATA BXY(1525),BXY(1526),BXY(1527),BXY(1528)/-64,  0,  4,  1/
      DATA BXY(1529),BXY(1530),BXY(1531),BXY(1532)/  5,  2,  5,  4/
      DATA BXY(1533),BXY(1534),BXY(1535),BXY(1536)/  3,  5,  3,  6/
      DATA BXY(1537),BXY(1538),BXY(1539),BXY(1540)/  4,  0,  0,  0/
      DATA BXY(1541),BXY(1542),BXY(1543),BXY(1544)/  8,  6,  8,  6/
      DATA BXY(1545),BXY(1546),BXY(1547),BXY(1548)/  0,  0,  0,  6/
      DATA BXY(1549),BXY(1550),BXY(1551),BXY(1552)/  8,-64,  0,  8/
      DATA BXY(1553),BXY(1554),BXY(1555),BXY(1556)/  6,  0,  0,  0/
      END

      block data colors
c------------------------------------------------------------
c
c     Purpose: this data is the font  for the hpgl colors
c
c------------------------------------------------------------
      USE hpgl1
      implicit double precision (a-h,o-z)
      common /psplot1/ scl,xoff,yoff

      data scl,xoff,yoff /400.,50.,50./

      end

