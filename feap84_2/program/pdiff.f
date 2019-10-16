c$Id:$
      function pdiff(x,i,ndm,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute the difference between maximum and minimum
c               nodal coordinates in direction-i.

c      Inputs:
c         x(ndm,* ) - Nodal coordinates for mesh
c         i         - Direction of comparison
c         ndm       - Spatial dimension of mesh
c         numnp     - Number of nodes in mesh

c      Outputs:
c         pdiff     - Difference between maximum and minimum
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   i, n, ndm, numnp
      real*8    pdiff,  xmx, xmn, x(ndm,numnp)

      save

      do n = 1,numnp
         if(mr(np(190)+n-1).ge.0) go to 110
      end do ! n

      if(ior.gt.0) then
        write(iow,3000)
        write(ilg,3000)
      else
        write(*,3000)
      endif
      pdiff = 0.0d0
      return

110   xmx = x(i,n)
      xmn = x(i,n)
      do n = 1,numnp
        if(mr(np(190)+n-1).ge.0) then
          xmx = max(xmx,x(i,n))
          xmn = min(xmn,x(i,n))
        endif
      end do ! n

      pdiff = xmx - xmn

c     Format

3000  format(' *ERROR* PDIFF: Coodinates are unspecified')

      end
