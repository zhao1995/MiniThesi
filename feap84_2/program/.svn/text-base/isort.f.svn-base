c$Id:$
      subroutine isort(ifeld,n)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Sorts the entries in array ifeld into increaseing order.
c               Required by sparse solver only.

c     Inputs:
c        ifeld(n) - Array of integers to sort
c        n        - Number of entries to sort

c     Output:
c        ifeld(n) - Sorted array of integers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   flag
      integer   i,ihilf,n
      integer   ifeld(n)

      save

c     Sort integer array using bubble sort algorithm

      flag = .true.
      do while( flag )
        flag = .false.
        do i = 1, n-1
          if (ifeld(i).gt.ifeld(i+1)) then
            flag       = .true.
            ihilf      = ifeld(i)
            ifeld(i)   = ifeld(i+1)
            ifeld(i+1) = ihilf
          end if
        end do ! i
      end do ! while flag = true

      end
