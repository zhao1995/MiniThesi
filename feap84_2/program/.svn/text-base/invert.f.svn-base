c$Id:$
      subroutine invert(a,nmax,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Invert small square matrix

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
      real*8    d, a(ndm,*)

      save

      do n = 1,nmax
        if(a(n,n).ne.0.0d0) then
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
c       else
c         write(*,*) ' *WARNING* Zero pivot in INVERT'
        endif
      end do ! n

      end
