c$Id:$
      subroutine fpplcl()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Write postscript filename to output file         29/12/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Close a PostScript file.

c      Inputs:
c         none

c      Outputs:
c         none      - Outputs written to postscript file: feappost.-
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pdatps.h'
      include  'plflag.h'
      include  'psdat4.h'
      include  'psdat6.h'

      character llx*9,lly*9,urx*9,ury*9
      integer   nc

      save

c     Close line with stroke if necessary

      if(lstrk) then
        call fppsin('s')
      endif

c     Add closing information to file

      call fppsdu()
      call fppsin('gr showpage')
      call fppsdu()
      call fppsin('%%Trailer')
      call fppsdu()
      call fppsin('%%EOF')
      call fppsdu()

c     Convert bounding box coordinates to character array

      if(psfram) then
        write(llx,'(i9)') -nint(yur*pscal) - 3 + 625
        write(lly,'(i9)')  nint(xll*pscal) - 3 -  10
        write(urx,'(i9)') -nint(yll*pscal) + 3 + 625
        write(ury,'(i9)')  nint(xur*pscal) + 3 -  10
      else
        write(llx,'(i9)')  nint(xll*pscal) - 3
        write(lly,'(i9)')  nint(yll*pscal) - 3
        write(urx,'(i9)')  nint(xur*pscal) + 3
        write(ury,'(i9)')  nint(yur*pscal) + 3
      endif

c     Create 'Feap#.eps file

      call feapbb(fname,llx,lly,urx,ury)

      nc = index(fname,' ')
      if(nc.eq.0) then
        nc = 17
      endif
      if(ior.lt.0) write(*,2000) fname(1:nc)
      write(iow,2001) fname(1:nc)

2000  format(' --> Closing FEAP PostScript file: ',a)
2001  format(' --> PostScript File: ',a)

      end
