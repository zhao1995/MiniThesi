c$Id:$
      subroutine chlfac(s,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Factor a positive definite matrix using Cholesky
c               Method.  Array stored by columns for upper part.

c      Inputs:
c         s(*) - Unfactored array
c         nn     - Size of array

c      Outputs:
c         s(*,*) - Factored array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,id,jd,nn
      real*8    s(*), dot

      save

c     Choleski factorization of a symmetric, positive definite matrix

      s(1) = 1.d0/sqrt(abs(s(1)))
      jd   = 1
      do j = 2,nn
        id = 0
        do i = 1,j-1
          if(i.gt.1) s(jd+i) = s(jd+i) - dot(s(id+1),s(jd+1),i-1)
          id = id + i
          s(jd+i) = s(jd+i)*s(id)
        end do ! i
        s(jd+j) = 1.d0/sqrt(abs(s(jd+j) - dot(s(jd+1),s(jd+1),j-1)))
        jd = jd + j
      end do ! j

      end
