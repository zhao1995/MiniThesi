c$Id:$
      subroutine formhh(hh,g,neqg,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Forms finite element G and H arrays

c      Inputs:
c         neqg   - Number of constraint equations
c         neq    - Number of active equations

c      Outputs:
c         hh(neqg,*) - Export matrix  (TEMP5: hr(np(115)))
c         g(neq,*,2) - G-vector       (TEMP6: hr(np(116)))
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'compas.h'
      include  'ndata.h'
      include  'part0.h'
      include  'pointer.h'
      include  'print.h'
      include  'p_int.h'
      include  'comblk.h'

      integer   neqg,neq,mm,nn, i
      real*8    hh(neqg,neqg),g(neq,neqg,2)

      save

c     Loop through equations to form Gu' = A_inv * Gu

      do mm = 1,neqg

c       Solve linear equations with Gu as right-hand side

        fp(1) = na
        fp(2) = nau
        fp(3) = nal
        fp(4) = np(20+npart)
        call psolve(ittyp,g(1,mm,1),fp,.false.,.true.,.true.,prnt)

        do nn = 1,neqg
          do i = 1,neq
            hh(nn,mm) = hh(nn,mm) - g(i,nn,2)*g(i,mm,1)
          end do ! i
        end do ! nn
      end do ! mm

      end
