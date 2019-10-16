c$Id:$
      subroutine propcg(jc,ir,jp,neq,ittyp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute profile for PCG preconditioners

c      Inputs:
c         jc(*)  - Pointer to ends of row/columns in compact store
c         ir(*)  - Location of non-zero components in each row/col
c         neq    - Number of equations
c         ittyp  - Maximum row/column length for preconditioner

c      Outputs:
c         jp(*)  - Pointers to ends of profile row/columns
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,i1,i2,ittyp, n, neq, jc(*), ir(*), jp(*)

      save

      do n = 1,neq
        jp(n) = 0
      end do ! n

      do n = 2,neq
        i1 = 0
        do i = jc(n-1)+1,jc(n)
          i2 = n - ir(i)
          if(i2.gt.0. and. i2.le.ittyp) then
            i1 = max(i1,i2)
          endif
        end do ! i
        jp(n) = max(jp(n),i1)
      end do ! n

c     Convert to profile pointers

      do i = 2,neq
        jp(i) = jp(i) + jp(i-1)
      end do ! i

      write(*,*) ' STORAGE FOR AUR =',jp(neq)

      end
