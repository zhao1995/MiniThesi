c$Id:$
      subroutine cdured(alr,ali,aur,aui,adr,adi,jh, djr,dji)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Reduce diagonal element in unsymmetric triangular
c               decomposition for complex matrix stored in real form

c      Inputs:
c         alr(*) - Real      part of lower terms
c         ali(*) - Imaginary part of lower terms
c         aur(*) - Real      part of upper terms
c         aui(*) - Imaginary part of upper terms
c         adr(*) - Real      part of diagonal terms
c         adi(*) - Imaginary part of diagonal terms
c         jh     - Length of row/column

c      Outputs:
c         djr    - Real      part of reduced diagonal
c         dji    - Imaginary part of reduced diagonal
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   j,jh
      real*8    udr,udi,djr,dji
      real*8    alr(jh),aur(jh),adr(jh), ali(jh),aui(jh),adi(jh)

      save

c     Computation for unsymmetric matrices
c     Scale upper U vector by reciprocal diagonals D ( forms U )

      do j = 1,jh
        udr    = aur(j)*adr(j) - aui(j)*adi(j)
        udi    = aur(j)*adi(j) + aui(j)*adr(j)
        aur(j) = udr
        aui(j) = udi
      end do ! j
c                        _
c     Dot product of L * U

      do j = 1,jh
        djr    = djr - aur(j)*alr(j) + aui(j)*ali(j)
        dji    = dji - aur(j)*ali(j) - aui(j)*alr(j)
      end do ! j
c                                                            _
c     Scale lower L vector by reciprocal diagonals D ( forms L )

      do j = 1,jh
        udr    = alr(j)*adr(j) - ali(j)*adi(j)
        udi    = alr(j)*adi(j) + ali(j)*adr(j)
        alr(j) = udr
        ali(j) = udi
      end do ! j

      end
