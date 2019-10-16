c$Id:$
      logical function ralloc(num,name,unit)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Read a dictionary entry for restart.
c               Pointer defined for integer (single) and real
c               (double precision arrays).

c      Inputs:
c         num        - Entry number for array (see below)
c         name       - Name of array          (see below)
c         unit       - Logical unit to use for a read

c      Outputs:
c         np(num)    - Pointer to first word of array memory
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'iofile.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    pcomp, setvar, palloc
      character  name*(*), dname*5
      integer    num, unit, length, precis, i

      save

c     Read length and precision

      read(unit) dname, length, precis

c     For arrays with non-zero length input the record

      if(length.gt.0) then
        if(pcomp(dname,name,5)) then

c         Allocate storage for array

          setvar = palloc(num,dname,length,precis)

c         Input array

          if(precis.eq.1) then
            read(unit) (mr(np(num)+i),i=0,length-1)
          elseif(precis.eq.2) then
            read(unit) (hr(np(num)+i),i=0,length-1)

c         Error on precision

          else
            write(iow,3000) num,name
            if(ior.lt.0) then
              write(*,3000) num,name
            endif
            call plstop()
          endif

c       Error on name

        else
          write(iow,3001) num, name,dname
          if(ior.lt.0) then
            write(*,3001) num, name,dname
          endif
          call plstop()
        endif
      endif

      ralloc = .true.

c     Format

3000  format(' **RESTART ERROR** Incorrect length or precision for',
     &       '                   Array Number',i4,' Named ',a/)

3001  format(' **RESTART ERROR** Incorrect name for array number',i4/
     &       '                   Request Name; ',a/
     &       '                   Input   Name; ',a/)

      end
