c$Id:$
      subroutine pformrx(id,ixt,ndf,numnp,neqms)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine number of flexible equations

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'part0.h'
      include  'part1.h'

      integer   n,ndf,numnp,neqms, j
      integer   id(ndf,*),ixt(*)

      save

      neqms = 0
      do n = 1,numnp
        if(ixt(n).eq.0) then
          do j = 1,ndf
            if(ndfp(j).eq.npart) then
              neqms = max(neqms,id(j,n))
            end if
          end do ! n
        end if
      end do ! n

      end
