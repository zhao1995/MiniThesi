c$Id:$
      subroutine memchk

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Limit memory use when maxuse > 0

c      Inputs:
c         none

c      Outputs:
c         none      - Program stops when memory use exceeded
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotd.h'
      include  'iofile.h'
      include  'memuse.h'

      character c*1
      integer   i

c     Output dictionary names

      totimem = 0
      totrmem = 0
      do i = 1,ndict
        if(iprec(i).eq.1) then
          totimem = totimem + ipoint(i)
        else
          totrmem = totrmem + ipoint(i)
        endif
      end do ! i
      if(totrmem + totimem.gt.maxuse) then
        write(  *,2000) totimem,totrmem,maxuse
        write(iow,2000) totimem,totrmem,maxuse
        write(ilg,2000) totimem,totrmem,maxuse
        read(*,'(a)') c
        call plstop()
      endif

c     Formats

2000  format(10x,'Total memory used by FEAP:'/
     &       20x,'Integer Arrays  = ',1i9/
     &       20x,'Real    Arrays  = ',1i9/
     &       20x,'Maximum allowed = ',1i9/ 10x,'Press Enter')

      end
