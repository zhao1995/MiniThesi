c$Id:$
      subroutine rfile(iunit, fname, irec, n, a)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Read A(n) array from file = 'fname',

c      Inputs:
c         iunit - Logical unit number of read
c         fname - File name for read
c         irec  - Block number for read
c         n     - Number of values to input

c      Outputs:
c         a(*)  - Block of matrix
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      include  'iofile.h'

      character fname*(*), tfile*12
      integer   iunit, irec, n,nm, i,j
      real*8    a(n)

      save

c     Open file for read - sequential access.

      if(irec.lt.10) then
        write(tfile,'(a6,i1)') fname,irec
      elseif(irec.lt.100) then
        write(tfile,'(a6,i2)') fname,irec
      else
        write(tfile,'(a5,i3)') fname,irec
      endif

      open(unit= iunit, file=tfile, access = 'sequential',
     &     form = 'unformatted', status = 'unknown',
     &     iostat = nm, err = 1000)

c     Read a from file.

      do i = 1,n,8000
        j = min(8000,n-i+1)
        call ioblok(iunit,a(i),j, 1)
      end do ! i

c     Close file.

      close(iunit)

      return

1000  write(*,2000)
      write(iow,2000)
2000  format(' *ERROR* Unexpected error in RFILE')

      end
