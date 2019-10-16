c$Id:$
      subroutine ploade(id,fn,u, work)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Adds work by external forces at current time

c      Inputs:
c         id(*)    - Equation numbers for each degree of freedom
c         fn(*)    - Force at t_n and t_n+1
c         u(*)     - Solution at t_n+1

c      Outputs:
c         work     - Work of external nodal forces on flexible body
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'part0.h'
      include  'sdata.h'

      integer   i,j,n, id(nneq)
      real*8    work,  fn(nneq),u(nneq)

      save

c     Compute work by external loads

      do i = 1,ndf ! {
        if(ndfp(i).eq.npart) then
          do n = i,nneq,ndf ! {
            j = id(n)
            if(j.gt.0) work = work + fn(n)*u(n)
          end do ! n          }
        endif
      end do ! i     }

      end
