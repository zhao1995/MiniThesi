c$Id:$
      subroutine sortpt(sv,npv)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sort list of points for dplot/splot display

c      Inputs:
c         sv(2,*)   - List of unsorted values
c         npv       - Number of values

c      Outputs:
c         sv(2,*)   - List of sorted values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,n,npv
      real*8    sv(2,npv),sn

      save

c     Sort on first values

      do n = 1,npv-1
        do i = n+1,npv
          if(sv(1,i).lt.sv(1,n)) then
            sn = sv(1,n)
            sv(1,n) = sv(1,i)
            sv(1,i) = sn
            sn = sv(2,n)
            sv(2,n) = sv(2,i)
            sv(2,i) = sn
          endif
        end do ! i
      end do ! n

      end
