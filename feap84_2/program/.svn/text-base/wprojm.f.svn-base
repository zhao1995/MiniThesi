c$Id:$
      subroutine wprojm(a,nn,ah)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Outputs projected subspace arrays: G and H

c      Inputs:
c         a(*)        - Array to output
c         nn          - Number row/columns in array
c         ah          - Header to write (G or H)

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      character ah*(*)
      integer   nn,i,k,n
      real*8    a(*)

      save

      i = 0
      write(iow,2000) ah
      do n = 1,nn
        write(iow,2001) (a(i+k),k=1,n)
        i = i + n
      end do ! n
      if(ior.lt.0) then
        i = 0
        write(  *,2000) ah
        do n = 1,nn
          write(  *,2001) (a(i+k),k=1,n)
          i = i + n
        end do ! n
      endif

c     Formats

2000  format(' ',a,'-Matrix ')

2001  format(1p,8d10.2)

      end
