c$Id:$
      subroutine ioblok(iunit,a,j, i)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Unformatted I/O for data block

c      Inputs:
c         iunit    - Logical unit number for I/O
c         a(j)     - Block of data to I/O
c         j        - Number of terms in a
c         i        - Read if 1; else write

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   iunit,i,j
      real*8    a(j)

      save

      if(i.eq.1) then
        read (iunit) a
      else
        write(iunit) a
      endif

      end
