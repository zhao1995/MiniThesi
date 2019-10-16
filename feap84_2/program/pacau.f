c$Id:$
      subroutine pacau(icp, neq, i, j, ac, ir, jc, jp, au)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Copy compact form into profile form.

c      Inputs:
c         icp        - size of compacted array.
c         neq        - number of active equations.
c         i          - number of first equation in block
c         j          - number of last  equation in block
c         ac(icp)    - compacted coefficient array
c         ir(icp)    - row number of nonzero term in compacted column.
c         jc(neq)    - end of entries in ir for column
c         jp(icp)    - end of columns in au/al.

c      Outputs:
c         au(icp)    - block to be stored on disk.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   icp,neq, i,j, np0,np,ni,nj, ir(icp), jc(neq), jp(neq)
      real*8    ac(icp), au(icp)

      save

c     Zero au

      do ni = 1,icp
        au(ni) = 0.0d0
      end do ! ni

c     Enter nonzero terms into au

      np0 = jp(i-1) - 1
      do ni = i,j
        np = jp(ni) - ni - np0
        do nj = jc(ni-1)+1,jc(ni)
          if(ac(nj).ne.0.0d0.and.np+ir(nj).gt.0
     &      .and.ir(nj).ne.0.and.np+ir(nj).le.icp) then
            au(np+ir(nj)) = au(np+ir(nj)) + ac(nj)
          endif
        end do ! nj
      end do ! ni

      end
