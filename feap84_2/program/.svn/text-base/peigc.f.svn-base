c$Id:$
      subroutine peigc(s,ld,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compress element matrix for repeated dof equations for
c                eigen extraction

c      Inputs:
c         s(nst,nst)   - Element array uncompressed
c         ld(nst)      - List of equation numbers for dof's
c         nst          - Size of element array

c      Outputs:
c         s(nst,nst)   - Element array compressed
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    nst, ld(*), i,j,k
      real*8     s(nst,nst)

      save

c     Compress matrix for eigen extraction

      do i = 1,nst-1
        do j = i+1,nst
          if(ld(j).eq.ld(i)) then
            do k = 1,nst
              s(k,i) = s(k,i) + s(k,j)
              s(k,j) = 0.0d0
            end do ! k
            do k = 1,nst
              s(i,k) = s(i,k) + s(j,k)
              s(j,k) = 0.0d0
            end do ! k
          endif
        end do ! j
      end do ! i

      end
