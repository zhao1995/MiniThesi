c$Id:$
      subroutine cdsred(aur,aui,adr,adi,jh, djr, dji)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Reduce diagonal element in symmetric triangular
c               decomposition for complex matrix stored in real form

c      Inputs:
c         aur(*) - Real      part of upper terms
c         aui(*) - Imaginary part of upper terms
c         adr(*) - Real      part of diagonal terms
c         adi(*) - Imaginary part of diagonal terms
c         jh     - height of columns

c      Outputs:
c         djr    - Real      part of reduced diagonal
c         dji    - Imaginary part of reduced diagonal
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   j,jh
      real*8    udr,udi, djr,dji, aur(jh),adr(jh), aui(jh),adi(jh)

      save

c     Computation for symmetric matrices ( forms U )

      do j = 1,jh
        udr    = aur(j)*adr(j) - aui(j)*adi(j)
        udi    = aur(j)*adi(j) + aui(j)*adr(j)
        djr    = djr - aur(j)*udr + aui(j)*udi
        dji    = dji - aur(j)*udi - aui(j)*udr
        aur(j) = udr
        aui(j) = udi
      end do ! j

      end
