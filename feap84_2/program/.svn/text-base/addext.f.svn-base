c$Id:$
      subroutine addext(fnam,fext,ifnam,ifext)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:       Add extender to file name for disk i/o

c      Inputs:
c        fname(*)   - File name without extender
c        fext(*)    - Extender
c        ifnam      - Length of 'fname'
c        ifext      - Length of 'fext'

c      Outputs:
c        fname(*)   - File name with added  '.' and 'fext'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   ipos, iposl, iposx, ifnam,ifext
      character fnam*(*),fext*(*)

      save

      iposl = ipos(fnam,ifnam) + 1
      iposx = ipos(fext,ifext)

      if(iposl+iposx.gt.ifnam) then
        write(  *,*) ' File names too long. Require: ',iposl+iposx
        write(iow,*) ' File names too long. Require: ',iposl+iposx
        call plstop()
      else
        fnam(iposl:ifnam)         = '. '
        fnam(iposl+1:iposl+iposx) = fext(1:iposx)
      endif

      end
