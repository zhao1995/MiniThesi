c$Id:$
      subroutine wfile(iunit, fname, irec, n, a)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Writes A(n) array to file = 'fname+irec'

c      Inputs:
c         iunit     - Logical unit number for write
c         fname     - Name of file to open
c         irec      - File number to write
c         n         - Length of record to write
c         a(*)      - Array to write to disk

c      Outputs:
c         none      - Outputs go to disk
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      character fname*(*), tfile*12
      integer   iunit, irec, n, nm, i,j
      real*8    a(n)

      save

c     Open file to write to - sequential access.

      if(irec.lt.10) then
        write(tfile,'(a6,i1)') fname,irec
      elseif(irec.lt.100) then
        write(tfile,'(a6,i2)') fname,irec
      else
        write(tfile,'(a5,i3)') fname,irec
      endif

      open(unit= iunit, file=tfile, access = 'sequential',
     +     form = 'unformatted', status = 'unknown',
     +     iostat = nm, err = 1000)

c     Write A to file.

      do i = 1,n,8000
        j = min(8000,n-i+1)
        call ioblok(iunit,a(i),j, 2)
      end do ! i

c     Close file.

      close(iunit)

      return

1000  write(  *,2000) tfile
      write(iow,2000) tfile
      write(ilg,2000) tfile

2000  format(' *ERROR* WFILE: Unexpected error in ',a)

      end
