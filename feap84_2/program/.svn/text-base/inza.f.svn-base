c$Id:$
      integer function inza(n1, n2, ir, k, n)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Locate term in sparse matrix

c      Inputs:
c         n1     - First location to search in ir
c         n2     - Last location to search in ir
c         ir(*)  - Array with locations in row
c         k      - Column number to find
c         n      - Equation number

c      Outputs:
c         inza   - Location of term
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   n1,n2, i,k,n, ir(*)

      save

c     Find term for assembly.

      do i = n1, n2
        if(ir(i).eq. k) then
          inza = i
          return
        endif
      end do ! i

c     ERROR if loop exits.

      inza = 0  ! Error
      write(  *,2000) n,k,n1,n2,(ir(i),i=n1,n2)
      write(iow,2000) n,k,n1,n2,(ir(i),i=n1,n2)
      call plstop()

c     Format

2000  format(/5x,'INZA: Eq. No. =',i8,' Col. No. =',i8/
     &        5x,'     Col. top =',i8,' Col. Bot =',i8/
     &        5x,'     Column entries'/ (5x,8i8))

      end
