c$Id:$
      subroutine elcntl(ie,ix,lagre, ic)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute array for performing assembly of stiffness

c      Inputs:
c          ie(*)    -  Element parameter entries
c          ix(*)    -  Element conectivity array.
c          lagre(*) -  Lagrange multiplier entries

c      Outputs:
c          ic(*)    -  Element lagrange multipler equations
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'sdata.h'

      integer    i,inn, n,nlm, ma
      integer    ie(nie,*), ix(nen1,*), lagre(*), ic(*)

      save

c     Lagrange multiplier elements

      do n = 1,numel
        if(ix(nen1-1,n).ge.0) then
          ma  = ix(nen1,n)       ! material set number
          nlm = ie(nie-8,ma)     ! gets number of element equations
          inn = lagre(ix(nen+4,n)) - 1
          do i = 1,nlm
            ic(inn+i) = ic(inn+i) + 1
          end do ! i
        endif
      end do ! n

      end
