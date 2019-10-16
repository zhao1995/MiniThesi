c$Id:$
      subroutine modify(p,s,dul,nsiz,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add modify of imaginary part by real b.c. value  18/08/2010
c       2. Remove check on ld and ld array                  17/03/2011
c       3. Add 'nsiz' to limit modify for some cases        15/03/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Modify element residual for effects of specified
c               boundary values.

c      Real:       p(i,1) = p(i,1) - s(i,j,1)*dul(j)
c      Imaginary:  p(i,2) = p(i,2) - s(i,j,2)*dul(j)

c      Inputs:
c         p(nst,*)     - Unmodified residual from element
c         s(nst,nst,*) - Element tangent array
c         dul(*)       - Value of specified solution increments
c         nst          - Dimension of element arrays

c      Outputs:
c         p(nst,*)     - Residual modified for effect of increments
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'complx.h'

      integer   nsiz,nst,i,j
      real*8    p(nst,*),s(nst,nst,*),dul(nst)

      save

c     Loop over columns

      do j = 1,nsiz

c       Loop over rows to modify equations

        do i = 1,nsiz
          p(i,1) = p(i,1) - s(i,j,1)*dul(j)       ! Real      modify
          if(cplxfl) then
            p(i,2) = p(i,2) - s(i,j,2)*dul(j)     ! Imaginary modify
          endif
        end do ! i

      end do ! j

      end
