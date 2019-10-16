c$Id:$
      subroutine setext(name,next, fext, flag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    20/12/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:       Add extender file number

c      Inputs:
c        name       - Name for extender
c        next       - Number of extender
c        flag       - Increment 'next' if .true.

c      Outputs:
c        fext(*)    - Extender with number 'next' added
c        next       - Next number
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'

      logical    flag
      integer    next
      character  fext*(*), name*(*)

c     Set extender name and initial number

      fext(1:4) =  name(1:4)
      fext(5:8) = '0000'
      if(next.le.9) then
        write(fext(8:8),'(i1)') next
      elseif(next.le.99) then
        write(fext(7:8),'(i2)') next
      elseif(next.le.999) then
        write(fext(6:8),'(i3)') next
      elseif(next.le.9999) then
        write(fext(5:8),'(i4)') next
      else
        write(ilg,3000) name
        write(iow,3000) name
        call plstop()
      endif
      if(flag) next = next + 1

c     Format

3000  format(' *ERROR* PMESH: More than 10000 files for ',a)

      end
