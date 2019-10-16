c$Id:$
      subroutine fppsdu()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Output string of characters to PostScript file

c      Inputs:
c        none       - Input through common /plpost/

c      Outputs:
c        none       - Outputs are written to PostScript file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iodata.h'
      include  'plpost.h'

      integer   i, first, last

      save

      if (nxtchr .gt. 0) then

c       Write to lun

        do first = 1,nxtchr
          if(buffer(first).ne.' ') go to 100
        end do ! first
        return
100     do last = nxtchr,first,-1
          if(buffer(last).ne.' ') go to 200
        end do ! last

200     write (lun,'(80a1)') (buffer(i), i=first,last)
        nxtchr = 0

c       Clear buffer

        do i = 1, 80
          buffer(i) = ' '
        end do ! i

      end if

      end
