c$Id:$
      subroutine cclred(aur,aui,xjr,xji,nn, br,bi)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Complex columnwise reduction for backsubstitution

c      Inputs:
c         aur(*) - Real part of column above diagonal
c         aui(*) - Imaginary part of column above diagonal
c         xjr    - Real part of solution term
c         xji    - Imaginary part of solution term
c         nn     - Length of column

c      Outputs:
c         br(*)  - Real part of reduced column
c         bi(*)  - Imaginary part of reduced column
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,nn
      real*8    xjr, xji, aur(*),br(*), aui(*),bi(*)

      save

      do n = 1,nn
        br(n) = br(n) - aur(n)*xjr + aui(n)*xji
        bi(n) = bi(n) - aur(n)*xji - aui(n)*xjr
      end do ! n

      end
