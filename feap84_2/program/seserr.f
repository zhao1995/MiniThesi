c$Id:$
      subroutine seserr(s,st,ix,nen,npstr,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Extract global projected quantity

c      Inputs:
c        st(*)   - Global quantity
c        ix(*)   - Element node numbers
c        nen     - Max nodes/element
c        npstr   - Number stress quantities
c        numnp   - Number global nodes

c      Outputs:
c        s(*)    - Element quantity
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    nen,npstr,numnp, n,nn,i, ix(*)
      real*8     s(nen,*),st(numnp,*)

      save

      do n = 1,nen
        nn = ix(n)
        if(nn.gt.0) then
          do i = 1,npstr
            s(n,i) = st(nn,i)
          end do ! i
        endif
      end do ! n

      end
