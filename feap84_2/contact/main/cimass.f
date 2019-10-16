c$Id:$
      subroutine cimass(lmass,id,lg,minv,eoffs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Set up inverse of local lumped mass array

c      Inputs:
c         lmass(*) - Lumped mass array
c         id(*)    - ID array
c         lg(*)    - Global-local node array

c      Outputs:
c         minv(*)  - Inverse of lumped mass in local system
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'sdata.h'

      real*8    lmass(*), minv(ndf,*)
      integer   id(ndf,*), lg(*), eoffs, n, i

c     Loop over master/slave nodal points

      do n=1,eoffs
        if(lg(n) .ne. 0) then
          do i=1,ndf
            if(id(i,lg(n)) .gt. 0) then
              minv(i,n) = 1.d0/lmass(id(i,lg(n)))
            else
              minv(i,n) = 0.d0
            endif
          end do ! i
        endif
      end do ! n

      end
