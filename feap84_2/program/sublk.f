c$Id:$
      subroutine sublk(n1,n2,n3,ix,ni,ne,nen1,ma)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate node connections for block of 3-d elements

c      Inputs:
c         n1,n2,n3  - Number of increments in each direction of block
c         ni        - Initial element number
c         ne        - Initial element number
c         nen1      - Dimension of ix array
c         ma        - Material set number for block

c      Outputs:
c         ix(nen1,*)- Block of element connections
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   me,ni,ne,nen1,ma,n1,n2,n3,n31,i,j,k,n
      integer   ix(nen1,*)

      save

c     Compute element connections

      if(ne.le.0) return
      me = ne - 1
      n31 = n3*n1
      do k = 1,n2-1
        do j = 1,n1-1
          n = n3*(j-1 + n1*(k-1)) + ni
          do i = 1,n3-1
            n = n + 1
            me = me + 1
            ix(nen1,me) = ma
            ix(1,me)    = n - 1
            ix(2,me)    = n
            ix(3,me)    = n + n3
            ix(4,me)    = n + n3  - 1
            ix(5,me)    = n + n31 - 1
            ix(6,me)    = n + n31
            ix(7,me)    = n + n31 + n3
            ix(8,me)    = n + n31 + n3 - 1
          end do ! i
        end do ! j
      end do ! k

      end
