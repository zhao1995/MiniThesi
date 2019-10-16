c$Id:$
      subroutine prlist( ilist, nlist )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Specify list of data for use in output by list.

c      Inputs:
c         none

c      Outputs:
c         ilist(*)  - List of numbers for selected outputs
c         nlist     - Number of items in ilist
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      character word*4
      logical   errck, pinput
      integer   ilist(100), nlist, n,nn
      real*8    td(8)

      save

      data      word/'node'/

c     Input list of items for printing

      nlist = 0
10    if(ior.lt.0) then
        write(*,3000)
        call pprint('  >')
      endif
      errck = pinput(td,8)
      if(errck) go to 10
      do n = 1,8
        nn = nint(td(n))
        if(nn.le.0) go to 200
        nlist = nlist + 1
        ilist(nlist) = nn
      end do ! n
      go to 10

c     Output list

200   write(iow,2000) (word, nn=1,min(8,nlist))
      if(ior.lt.0) then
        write(*,2000) (word, nn=1,min(8,nlist))
      endif
      do n = 1,nlist,8
        write(iow,2001) (ilist(nn),nn=n,min(nlist,n+7))
        if(ior.lt.0) then
          write(*,2001) (ilist(nn),nn=n,min(nlist,n+7))
        endif
      end do ! n

c     Formats

2000  format(/'     L i s t    o f    P r i n t   I t e m s'/
     &   2x,8(4x,a4:))

2001  format(2x,8i8)

3000  format(' Input: List of nodes (8/record)')

      end
