c$Id:$
      subroutine feapbb(filer,llx,lly,urx,ury)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Puts bounding box at beginning of file

c      Inputs:
c         llx        - Bounding box lower left  x-coord
c         lly        - Bounding box lower left  y-coord
c         urx        - Bounding box upper right x-coord
c         ury        - Bounding box upper right y-coord

c      Outputs:
c         filer      - Filename with postscript data and bounding box
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      include  'iodata.h'

      character llx*9,lly*9,urx*9,ury*9
      character line*80, boundbox*50, filer*(*)
      logical   eofile
      integer   ii

      save

c     Set up bounding box record

      boundbox( 1:14) ='%%BoundingBox:'
      boundbox(15:23) = llx
      boundbox(24:32) = lly
      boundbox(33:41) = urx
      boundbox(42:50) = ury

c     Open and rewind write file

      open(unit=ios,file=filer,status='unknown')
      rewind(lun)
      rewind(ios)

c     Read records from 'temp.eps' copy to 'Feap#.eps'

      eofile = .true.
      do while (eofile)
        read(lun,'(a)',end=200) line

c       Non bounding box records

        if(line(3:7).ne.'Bound') then
          do ii = 80,1,-1
            if(line(ii:ii).ne.' ') go to 100
          end do ! ii
100       write(ios,'(a)') line(1:ii)

c       BoundingBox record

        else
          write(ios,'(a50)') boundbox
        endif
      end do ! while

200   close(lun,status='delete')
      close(ios)

      end
