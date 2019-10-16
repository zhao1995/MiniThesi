c$Id:$
      subroutine pushr4(tl,tr,dm,ds,detf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Push forward 4th rank tensor
c         ds(i,j) = tl(k,i)*dm(k,l)*tr(l,j)/detf

c     Inputs:
c         dm(6,6) - material moduli
c         tl(6,6) - left  transformation array
c         tr(6,6) - right transformation array
c         detf    - determinant of deformation gradient
c     Outputs:
c         ds(6,6) - spatial moduli
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i,j,k
      real*8   detf, jrec, tl(6,6), tr(6,6), dm(6,6), ds(6,6), dtj(6)

c     Reciprocal deformation gradient determinant

      jrec = 1.d0/detf

c     Compute matrix product: dtj = dm*tr

      do j = 1,6
        do i = 1,6
          dtj(i) = 0.0d0
          do k = 1,6
            dtj(i) = dtj(i) + dm(i,k)*tr(k,j)
          end do ! k
        end do ! i

c       Compute spatial tensor: ds = tl_trans*dt

        do i = 1,6
          ds(i,j) = 0.0d0
          do k = 1,6
            ds(i,j) = ds(i,j) + tl(k,i)*dtj(k)
          end do ! k
          ds(i,j) = ds(i,j)*jrec
        end do ! i
      end do ! j

      end
