c$Id:$
      logical function walloc(num,name,length,precis,unit)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Write a dictionary entry for restart.
c               Pointer defined for integer (single) and real
c               (double precision arrays).

c      Inputs:
c         num        - Entry number for array (see below)
c         name       - Name of array          (see below)
c         length     - Length of array to save
c         precis     - Precision of array to save: 1 = integer; 2 = real
c         unit       - Logical unit to use for a write

c      Outputs:
c         none       - Information saved using 'unit'
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'iofile.h'
      include   'pointer.h'
      include   'comblk.h'

      character  name*(*),dname*5
      integer    num, unit, length, precis, i

      save

c     Read length and precision

      dname = name
      write(unit) dname, length, precis

      if(precis.eq.1) then
        write(unit) (mr(np(num)+i),i=0,length-1)
      elseif(precis.eq.2) then
        write(unit) (hr(np(num)+i),i=0,length-1)
      else
        write(iow,3000) num,name
        if(ior.lt.0) then
          write(*,3000) num,name
        endif
        call plstop()
      endif

      walloc = .true.

c     Format

3000  format(' **RESTART ERROR** Incorrect length or precision for',
     &       '                   Array Number',i4,' Named ',a/)

      end
