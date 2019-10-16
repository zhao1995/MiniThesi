c$Id:$
      subroutine fppsop(scal)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Open a new PostScript file to receive plot data
c               Maximum files: 456,976 (FeapAAAA.eps to FeapZZZ.eps)

c      Inputs:
c         scal      - Scale factor for plot data to be written

c      Outputs:
c         none      - Outputs are written to PostScript file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'pdatps.h'
      include  'pfeapb.h'
      include  'plflag.h'
      include  'plpost.h'
      include  'psdat2.h'
      include  'psdat4.h'
      include  'psdat5.h'
      include  'psdat6.h'

      logical   fexist
      integer   i, ii, iln(2), nch
      real*8    scal
      character cdate*24, uname*8, string*20, title*68, ftemp*17

      save

      nxtchr = 0

c     Get current date and user's name for banner page

      call fdate(cdate)
      call getlog(uname)

c     Find a file name that does not already exist

      if(pfeap_on) then
        call fparsop(fname,ftemp,nch)
      else
        fname = 'FeapAAAA.eps'
        ftemp = 'temp.eps'
        nch   =  8
      endif
      inquire(file = fname,exist = fexist)
      do while(fexist)
        call postname(fname,nch)
        inquire(file = fname,exist = fexist)
      end do ! while

c     Set initial BoundingBox coordinates

      xll = 5800
      yll = 4500
      xur = 0
      yur = 0

c     Add problem title to file

      ii = 0
      do i = 1,67,4
        ii           = ii + 1
        title(i:i+3) = head(ii)
      end do ! i

c     Show initialization worked, i.e. we opened file.

      if(ior.lt.0) write(*,2000) fname(1:nch+4)

      open(unit=lun,file=ftemp,status='unknown')

c     Write header information to file

      call fppsin('%!PS-Adobe-3.0 EPSF-3.0')
      call fppsdu()
      call fppsin('%%BoundingBox: (atend) ')
      call fppsdu()
      call fppsin('%%Title: '//title)
      call fppsdu()
      call fppsin('%%Creator: '//uname)
      call fppsdu()
      call fppsin('%%Creation Date: '//cdate//' ')
      call fppsdu()
      call fppsin('%%EndComments ')
      call fppsdu()

c     Set procedure definitions

      call fppsin('/m {moveto} bind def /l {lineto} bind def ')
      call fppsdu()
      call fppsin('/s {stroke} bind def /f {fill} bind def ')
      call fppsdu()
      call fppsin('/n {newpath} bind def /c {closepath} bind def ')
      call fppsdu()
      call fppsin('/g {setgray} bind def /h {setrgbcolor} bind def ')
      call fppsdu()
      call fppsin('/d {setdash} bind def /lw {setlinewidth} bind def ')
      call fppsdu()

c     Set gray scale shades

      call fppsin('/G0 { 0.0 g} bind def /G1 { 1.0 g} bind def ')
      call fppsdu()
      call fppsin('/G2 {0.95 g} bind def /G3 {0.88 g} bind def ')
      call fppsdu()
      call fppsin('/G4 {0.81 g} bind def /G5 {0.74 g} bind def ')
      call fppsdu()
      call fppsin('/G6 {0.67 g} bind def /G7 {0.60 g} bind def ')
      call fppsdu()
      call fppsin('/G8 {0.53 g} bind def /G9 {0.46 g} bind def ')
      call fppsdu()
      call fppsin('/Ga {0.39 g} bind def /Gb {0.32 g} bind def ')
      call fppsdu()
      call fppsin('/Gc {0.25 g} bind def /Gd {0.18 g} bind def ')
      call fppsdu()
      call fppsin('/Ge {0.11 g} bind def /Gf {0.04 g} bind def ')
      call fppsdu()
      call fppsin('/Gg {0.01 g} bind def /Gh {0.00 g} bind def ')
      call fppsdu()

c     Set color scale hues

      call fppsin('/h0 { 0.0 0.0 0.0 h} bind def') ! black
      call fppsdu()
      call fppsin('/h1 { 1.0 1.0 1.0 h} bind def') ! white
      call fppsdu()
      call fppsin('/h2 { 1.0 0.0 0.0 h} bind def') ! red
      call fppsdu()
      call fppsin('/h3 { 0.0 1.0 0.0 h} bind def') ! green
      call fppsdu()
      call fppsin('/h4 { 0.0 0.0 1.0 h} bind def') ! blue
      call fppsdu()
      call fppsin('/h5 { 1.0 1.0 0.0 h} bind def') ! yellow
      call fppsdu()
      call fppsin('/h6 { 0.0 1.0 1.0 h} bind def') ! cyan
      call fppsdu()
      call fppsin('/h7 { 1.0 0.0 1.0 h} bind def') ! magenta
      call fppsdu()
      call fppsin('/h8 { 1.0 .65 0.0 h} bind def') ! orange
      call fppsdu()
      call fppsin('/h9 { 1.0 0.5 .31 h} bind def') ! coral
      call fppsdu()
      call fppsin('/ha { .67 1.0 .18 h} bind def') ! green yellow
      call fppsdu()
      call fppsin('/hb { .96 .87 .70 h} bind def') ! wheat
      call fppsdu()
      call fppsin('/hc { .25 .41 1.0 h} bind def') ! royal blue
      call fppsdu()
      call fppsin('/hd { .63 .125 .94 h} bind def')! purple
      call fppsdu()
      call fppsin('/he { 0.5 1.0 .83 h} bind def') ! aquamarine
      call fppsdu()
      call fppsin('/hf { .75 .75 .75 h} bind def') ! gray
      call fppsdu()
      call fppsin('/hg { .83 .83 .83 h} bind def') ! light gray
      call fppsdu()
      call fppsin('/hh { 0.0 0.0 0.0 h} bind def') ! black
      call fppsdu()

c     Set Line types

      call fppsin('/l1 { [] 0 d } bind def ')
      call fppsdu()
      call fppsin('/l2 { [5 30] 0 d } bind def ')
      call fppsdu()
      call fppsin('/l3 { [40 20 5 20] 0 d } bind def ')
      call fppsdu()
      call fppsin('/l4 { [40] 0 d } bind def ')
      call fppsdu()
      call fppsin('/l5 { [60] 0 d } bind def ')
      call fppsdu()
      call fppsin('/l6 { [5 20 5 40 40 40] 0 d } bind def ')
      call fppsdu()
      call fppsin('/l7 { [40 60 80 60] 0 d } bind def ')
      call fppsdu()
      call fppsin('/l8 { [80] 0 d } bind def ')
      call fppsdu()

c     Set for 12 point type in landscape mode (10 in portrait)

      if(psfram) then
        call fppsin
     &   ('/H1 {/Helvetica findfont  91 scalefont setfont} bind def')
        call fppsin
     &   ('/H2 {/Helvetica findfont 146 scalefont setfont} bind def')
        call fppsin
     &   ('/H3 {/Helvetica findfont 218 scalefont setfont} bind def')
      else
        call fppsin
     &   ('/H1 {/Helvetica findfont 100 scalefont setfont} bind def')
        call fppsin
     &   ('/H2 {/Helvetica findfont 160 scalefont setfont} bind def')
        call fppsin
     &   ('/H3 {/Helvetica findfont 240 scalefont setfont} bind def')
      endif
      call fppsdu()
      call fppsin('/w {stringwidth pop 2 div neg 0 rmoveto} bind def')
      call fppsdu()

c     Set clipping definitions

      call fppsin('/gr {grestore} bind def')
      call fppsdu()
      call fppsin('/gs {gsave} bind def')
      call fppsdu()
      call fppsin('/cl {gr gs 802 802 3333 3333 rectclip} bind def')
      call fppsdu()
      call fppsin('/fl {gr gs   0   0 5800 4800 rectclip} bind def')
      call fppsdu()

c     End of prolog

      call fppsin('%%EndProlog ')
      call fppsdu()

c     Start landscape mode plot

      if(psfram) then
        call fppsin('%Landscape mode ')
        call fppsdu()
        if(blk) then
          call fppsin(' 0 0 0 h ')
          call fppsin('n 0 0 m 612 0 l 612 792 l 0 792 l c f')
          call fppsdu()
        end if
        call fppsin('90 rotate -10 -625 translate ')
        write(string,'(f6.4,f7.4,a7)') scal*0.1333,scal*0.1333,' scale '
        call fppsin(string)
        pscal = scal*0.1333

c     Start portrait mode plot

      else
        call fppsin('%Portrait mode ')
        call fppsdu()
        write(string,'(f6.4,f7.4,a7)') scal*0.1,scal*0.1,' scale '
        pscal = scal*0.1
        call fppsin(string)
        if(blk) then
          call fppsdu()
          call fppsin('0 0 0 h n ')
          call fppsin('318 318 m 5800 318 l 5800 4495 l 318 4495 l c f')
        end if
      endif
      call fppsdu()
      call fppsin('1 setlinecap 1 setlinejoin ')

c     Initialize plot state for lines, fills, and colors to false

      lstrk = .false.
      lfill = .false.
      xold  = -9980
      yold  = -9980
      dold  = -1
      lwold = -1
      clin  = 'G0'
      oclin = '  '
      cvar  = ' G0'
      ocvar = ' '
      colv  = 'z '
      ocolv = '  '

      iln(1) = 0
      iln(2) = 1
      call plline( iln )
      call fppsin(' gs n')
      call fppsdu()

c     Format

2000  format(' --> Opening FEAP PostScript file: ',a )

      end
