c$Id:$
      subroutine pgetd( name, point, lengt, prec , flag )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Retrieve an array from dictionary

c      Inputs:
c         name     - Name of array to retrieve

c      Outputs:
c         point    - Pointer to array
c         lengt    - Length of array
c         prec     - Precision of array
c         flag     - Flag, true if array found
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotn.h'
      include  'allotd.h'
      include  'debugs.h'
      include  'iofile.h'
      include  'pointer.h'

      include  'p_point.h'

      character name*(*),dname*5
      logical   pcomp,flag
      integer   lengt, prec, i

      save

c     Search dictionary for name

      dname = name

      do i = 1,ndict

c       Assign pointer, length, and precision

        if( pcomp(dname, dict(i), 5) ) then
          point =  np(dlist(i))
          lengt =  ipoint(i)
          prec  =  iprec(i)
          flag  = .true.
          return
        endif
      end do ! i

c     Output error message if not found

      if(debug) then
        if(ior.lt.0) write(*,2000) dname(1:5)
        write(iow,2000) dname(1:5)
c       if(ior.gt.0) call plstop()
      end if
      flag  = .false.
      point = 0
      lengt = 0
      prec  = 0

c     Format

2000  format(' *WARNING* Check for ',a5,
     &       ': Array not allocated for this problem.')

      end
