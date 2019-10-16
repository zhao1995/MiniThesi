c$Id:$
      subroutine udebug(string,iopt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  User solver interface

c     Inputs:
c       string      - Character string to output
c       iopt        - Integer

c     Outputs:
c       Write output to unit 'iunit'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'setups.h'

      logical    lopen
      character  string*(*), text*30, filnam*8
      integer    iopt,iunit

c     Set unit number

      iunit = 90 + rank

c     Check for existing output file

      inquire(unit = iunit, opened = lopen)

      if(.not.lopen) then
        filnam = 'Udebug0'
        if(rank.le.10) then
          write(filnam(8:8),'(i1)') rank
        else
          write(filnam(7:8),'(i2)') rank
        endif
        open(unit = iunit, file = filnam)
      endif

      text = '  '
      text = string

      write(iunit,2000) text,iopt
      call pflush(iunit)

c     Output format

2000  format(a,'I =',i5)

      end
