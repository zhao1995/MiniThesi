c$Id:$
      subroutine mkface(iblend,lblend)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form faces for 3-d Block

c      Inputs:
c        iblend(*)    - Blending function descriptor array

c      Outputs:
c        lblend(20,*) - Global blend array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i,j, iblend(*),lblend(20,*), lface(4,6)

      save

      data  lface/  1, 2, 3, 4,   5, 6, 7, 8,   1, 4, 8, 5,
     &              2, 3, 7, 6,   1, 5, 6, 2,   4, 8, 7, 3/

c     Assign generation increments for faces

      lblend(1,1) = iblend(1)
      lblend(2,1) = iblend(2)
      lblend(3,1) = 1

      lblend(1,2) = iblend(1)
      lblend(2,2) = iblend(2)
      lblend(3,2) = 1

      lblend(1,3) = iblend(2)
      lblend(2,3) = iblend(3)
      lblend(3,3) = 1

      lblend(1,4) = iblend(2)
      lblend(2,4) = iblend(3)
      lblend(3,4) = 1

      lblend(1,5) = iblend(3)
      lblend(2,5) = iblend(1)
      lblend(3,5) = 1

      lblend(1,6) = iblend(3)
      lblend(2,6) = iblend(1)
      lblend(3,6) = 1

c     Assign side numbers to faces

      do i = 1,6
        do j = 1,4
          lblend(j+10,i) = iblend(10+lface(j,i))
        end do ! j
      end do ! i

      end
