c$Id:$
      subroutine cmprot(mo,lam, numnp, lmr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      1. Add output of rotation lambda array               05/06/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set rotational array for finite rotation transformation

c      Inputs:
c        m0(*)     - Indicator for active rotation at node
c        lam(:,:)  - Uncompressed rotation array for node
c        numnp     - Number of nodes

c      Outputs:
c        lam(:,:)  - Compressed rotation array for node
c        lmr       - Number active nodal rotations
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include  'iofile.h'
      include  'print.h'

      integer  numnp, lmr, n,i
      integer  mo(numnp,2)
      real*8   lam(54,numnp)

c     Compress and output nodal lambda arrays

      if(prt) then
        write(iow,2000)
      endif
      lmr = 0
      do n = 1,numnp
        if(mo(n,1).ne.0) then
          lmr = lmr + 1
          mo(n,2) = lmr
          if(lmr.ne.n) then
            do i = 1,54
              lam(i,lmr) = lam(i,n)
            end do ! i
          endif
          if(prt) then
            write(iow,2001) n,lam(1,n),lam(4,n),lam(7,n),
     &                        lam(2,n),lam(5,n),lam(8,n),
     &                        lam(3,n),lam(6,n),lam(9,n)
          endif
        endif
      end do ! n

c     Formats

2000  format(/5x,'N o d a l   L a m b d a   A r r a y'//
     &       9x,'Node   v(1)        v(2)        v(3)')

2001  format(/i12,1p3e12.4/(12x,1p3e12.4))

      end
