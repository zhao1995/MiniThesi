c$Id:$
      subroutine pexpd(a,t,id,ndf,nneq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Expand compressed vector to include all values
c               removed by boundary conditions.

c      Inputs:
c         t(*)      - Compressed vector
c         id(ndf,*) - Equation number array
c         ndf       - Number dof/node
c         nneq      - Number in uncompressed vector

c      Outputs:
c         a(*)      - Uncompressed vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'part0.h'

      integer   i,n,ndf,nneq, id(*)
      real*8    a(*),t(*)

      save

c     Exand compressed temporary vector

      do i = 1,ndf
        if(npart.eq.ndfp(i)) then
          do n = i,nneq,ndf
            if(id(n).gt.0) then
              a(n) = t(id(n))
            else
              a(n) = 0.0d0
            endif
          end do ! n
        endif
      end do ! i

      end
