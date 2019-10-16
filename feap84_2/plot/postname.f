c$Id:$
      subroutine postname( name, nch )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set filename for postscript output
c               Maximum files: 456,976 (FeapAAAA.eps to FeapZZZ.eps)

c      Inputs:
c         name      - Default name
c         nch       - Length of filename

c      Outputs:
c         name      - New filename
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      character  name*(*)

      integer    n, nc, nch
      logical    add

      save

c     Initialize

      add = .true.

c     Check names

      n = nch
      do while( add .and. n.gt.1 )

c       Get a character from 'name'

        nc = ichar(name(n:n))

c       Check that it is less than a 'Z'

        if(nc.lt.90) then
          name(n:n) =  char(nc+1)
          add       = .false.

c       It is a 'Z' (or something erroneous!).  Do next column

        else
          name(n:n) = 'A'
        endif

        n = n - 1

      end do ! while

c     Too many files exist!

      if(n.le.1) then
        write(*,*) ' *ERROR* POSTNAM: Too many file names'
      endif

      end
