c$Id:$
      subroutine pwopn ()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Replace MSFLIB by DFLIB                          08/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      use DFLIB

      implicit none

      type(windowconfig)  :: textscreen
      logical             :: status
      integer             :: numxpix,numypix,numcols,numrows
      integer             :: numclrs,fontsz

      save

      open(unit = 0, file = 'CON')

c     Get Current Input Data Window Size

      status  = getwindowconfig(textscreen)
      numxpix = textscreen.numxpixels
      numypix = textscreen.numypixels
      numcols = textscreen.numtextcols
      numrows = textscreen.numtextrows
      numclrs = textscreen.numcolors
      fontsz  = textscreen.fontsize

c     Set up Input Data Window

      textscreen.numxpixels   =  640
      textscreen.numypixels   =  0.80*float(numypix)
      textscreen.numtextcols  =  80
      textscreen.numtextrows  =  36
      textscreen.numcolors    = -1
      textscreen.fontsize     = #0008000E
      textscreen.title = "F E A P     D a t a     I n p u t s"C
c     textscreen.bitsperpixel = -1

      status = setwindowconfig(textscreen)
      if(.not.status) status = setwindowconfig(textscreen)

      status = displaycursor($gcursoron)

      end
