c$Id:$
      subroutine reshis(ix,nen1,numel,n1, n2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change loop 'nh' limit to int4 difference        17/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Initialize t_n+1 history variables from final value
c               of variables at t_n

c      Inputs:
c         ix(nen1,*)  - Element connection/history pointer array
c         nen1        - Dimension of ix array
c         numel       - Number of elements in mesh
c         n1          - Pointer in ix to t_n   data
c         n2          - Pointer in ix to t_n+1 data

c      Outputs:
c         none        - Output is retained in pointers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'hdata.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   nen1,numel,n1,n2, n,nh, ix(nen1,numel)

      save

c     Move history variables from location 'n1' to 'n2'

      do n = 1,numel
        nh1 = np(49) + ix(n1,n) - 1
        nh2 = np(49) + ix(n2,n) - 1
        if(nh2.ne.nh1) then
          do nh = 1,abs(ix(n2,n) - ix(n1,n))
            hr(nh+nh2) = hr(nh+nh1)
          end do ! nh
        endif
      end do ! n

c     Move interface variables

      if(np(210).ne.0) then
        call ireshis(mr(np(210)), n1,n2)
      endif

      end
