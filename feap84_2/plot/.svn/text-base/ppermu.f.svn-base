c$Id:$
      subroutine ppermu(sym)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute reflection permutation values for graphics
c               reflections

c      Inputs:
c         sym(3)   - Symmetry indicator array

c      Outputs:
c         none     - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdatay.h'
      include  'pdatas.h'

      integer   i , j , k, sym(3)

      save

c     Initialize arrays

      call pzero(xsym,24)
      do i = 1,8
        isym(i) = i
        nfac(i) = nfac(1)
      end do ! i

c     Mark symmetries from inputs

      nsym = 0
      do i = 0,sym(1)
        do j = 0,sym(2)
          do k = 0,sym(3)
            if(iquad(i+1,j+1,k+1).eq.0) then
              nsym         = nsym + 1
              xsym(1,nsym) = 1.d0 - 2.d0*i
              xsym(2,nsym) = 1.d0 - 2.d0*j
              xsym(3,nsym) = 1.d0 - 2.d0*k
            end if
          end do ! k
        end do ! j
      end do ! i

c     Mark reflection properties (vsym=+1 for rhs system
c                                     =-1 for lhs system)
      do i = 1,nsym
        vsym(i) = 1.d0
        do j = 1,3
         vsym(i) = vsym(i)*xsym(j,i)
        end do ! j
      end do ! i

      lsym = 1

      end
