c$Id:$
      subroutine psetip(ip,numel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set integer array to its position value for sorting
c               elements by z-sort for hidden surface plots

c      Inputs:
c         numel    - Number of elements in mesh

c      Outputs:
c         ip(8,*)  - Integer array with positions set
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, n, numel, ip(8,numel)

      save

      do n = 1,numel
        do i = 1,8
         ip(i,n) = n
        end do ! i
      end do ! n

      end
