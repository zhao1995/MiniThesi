c$Id:$
      subroutine  sortjc(ir,jc,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Redefines row pointers and sorts column entries in
c              rows to ascending order
c     Inputs:
c         ir(neq)   - Row pointers from compro
c         jc(*)     - Unsorted column entries in rows
c         neq       - Number of equations
c     Outputs:
c         ir(neq+1) - Row pointers for sparse solver
c         jc(*)     - Sorted column entries
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i,neq
      integer  ir(*),jc(*)

      save

c     Expand row counter for sparse solver

      ir(neq+1) = 1
      do i = 1,neq
        ir(neq+i+1) = ir(i) + 1
      end do ! i

      do i = 1,neq
        call isort(jc(ir(neq+i)),ir(neq+i+1)-ir(neq+i))
      end do ! i

      end
