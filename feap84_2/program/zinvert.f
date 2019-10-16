c$Id:$
      subroutine zinvert(a,nmax,ndm,t)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Invert small square matrix, test for singularity

c      Inputs:
c         a(ndm,*) - Matrix to be inverted
c         nmax     - Size of upper submatrix to invert
c         ndm      - Dimension of array

c      Outputs:
c         a(ndm,*) - Submatrix replaces original terms, others not
c                    changed
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,n,ndm,nmax
      real*8    d, a(ndm,*), zpiv, tol, t(*)

      save

      data      tol /1.d-6/

      zpiv = a(1,1)*tol
      do n = 1,nmax
        t(n) = a(n,n)*tol
      end do ! n
      do n = 1,nmax
c       if(abs(a(n,n)).gt.zpiv) then
        if(abs(a(n,n)).gt.t(n)) then
          d = 1.d0/a(n,n)
          do j = 1,nmax
            a(n,j) = -a(n,j)*d
          end do ! j

          do i = 1,nmax
            if(n.ne.i) then
              do j = 1,nmax
                if(n.ne.j) a(i,j) = a(i,j) + a(i,n)*a(n,j)
              end do ! j
            endif
            a(i,n) = a(i,n)*d
          end do ! i
          a(n,n) = d
        else
          do i = 1,nmax
            a(i,n) = 0.0d0
            a(n,i) = 0.0d0
          end do ! i
        endif
      end do ! n

      end
