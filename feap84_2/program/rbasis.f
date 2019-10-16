c$Id:$
      subroutine rbasis(x, ebig)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set basis vectors for joint axes in reference state

c      Inputs:
c         ndm       - Spatial dimension of mesh
c         x(3,2)    - Coordinates defining direction of E_3

c      Outputs:
c         ebig(3,3) - Basis vectors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   i,imn
      real*8    enorm,tol, x(3,2),ebig(3,3)

      save

      data      tol /1.d-10/

c     Set initial vectors to zero

      do i = 1,3
        ebig(i,1) = 0.0d0
        ebig(i,2) = 0.0d0
        ebig(i,3) = 0.0d0
      end do ! i

c     Set axial vector: E_big(i,3)

      enorm = 0.0d0
      do i = 1,3
        ebig(i,3) = x(i,2) - x(i,1)
        enorm     = enorm + ebig(i,3)**2
      end do ! i
      if(enorm.gt.tol) then
        enorm = 1.d0/sqrt(enorm)
      else
        write(iow,3000) x
        call plstop()
      endif
      do i = 1,3
        ebig(i,3) = ebig(i,3)*enorm
      end do ! i

c     Find Minimum direction for orthogonal construction

      imn = 1
      do i = 2,3
        if(abs(ebig(i,3)).lt.abs(ebig(imn,3))) then
          imn = i
        end if
      end do ! i

c     Set unit vector in direction of minimum component

      ebig(imn,1) = 1.d0

c     Construct E_big(i,2) vector orthogonal to 3-vector

      ebig(1,2) = ebig(2,3)*ebig(3,1) - ebig(3,3)*ebig(2,1)
      ebig(2,2) = ebig(3,3)*ebig(1,1) - ebig(1,3)*ebig(3,1)
      ebig(3,2) = ebig(1,3)*ebig(2,1) - ebig(2,3)*ebig(1,1)
      enorm     = ebig(1,2)**2 + ebig(2,2)**2 + ebig(3,2)**2
      enorm     = 1.d0/sqrt(enorm)
      do i = 1,3
        ebig(i,2) = ebig(i,2)*enorm
      end do ! i

c     Construct E_big(i,1) vector orthogonal to 2 & 3-vectors

      ebig(1,1) = ebig(2,2)*ebig(3,3) - ebig(3,2)*ebig(2,3)
      ebig(2,1) = ebig(3,2)*ebig(1,3) - ebig(1,2)*ebig(3,3)
      ebig(3,1) = ebig(1,2)*ebig(2,3) - ebig(2,2)*ebig(1,3)
      enorm     = ebig(1,1)**2 + ebig(2,1)**2 + ebig(3,1)**2
      enorm     = 1.d0/sqrt(enorm)
      do i = 1,3
        ebig(i,1) = ebig(i,1)*enorm
      end do ! i

      call mprint(ebig,3,3,3,' E-Big')

3000  format(' *ERROR* Joint axis has zero length:'/
     &       5x,'X-1 = ',1p,3e12.4/
     &       5x,'X-2 = ',1p,3e12.4)
      end
