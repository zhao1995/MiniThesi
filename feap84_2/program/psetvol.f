c$Id:$
      subroutine psetvol(x,ndm,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add computation of center: 'xc0(3)'              13/03/2011
c       2. Remove volm0 and c0 from argument of routine     13/03/2011
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Compute volume of RVE

c     Inputs
c       x(ndm,numnp) - Nodal coordinates
c       nmd          - Problem space dimension
c       numnp        - Number of nodes on RVE

c     Output:          in common /elpers.h/
c       volm0        - Volume of RVE
c       xc(3)        - Side lengths coordinate
c       xc0(3)       - Mean coordinate
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'debugs.h'
      include   'elpers.h'
      include   'iofile.h'

      integer    ndm,numnp, n,i
      real*8     x(ndm,numnp), xmn,xmx

      save

      volm0 = 1.0d0
      do i = 1,ndm
        xmn = x(i,1)
        xmx = x(i,1)
        do n = 2,numnp
          xmn = min(xmn,x(i,n))
          xmx = max(xmx,x(i,n))
        end do ! n
        volm0  = volm0*(xmx - xmn)
        xc(i)  = (xmx - xmn)*0.5d0
        xc0(i) = (xmx + xmn)*0.5d0
      end do ! i

c     Output values of RVE

      write(iow,2000) volm0,(i,xc (i),i,xc0(i),i=1,ndm)
      if(debug) then
        write(*,2000) volm0,(i,xc (i),i,xc0(i),i=1,ndm)
      endif

2000  format(5x,'RVE Size Data'/
     &       10x,'Volume   =',1p,1e11.3/
     &      (10x,'Side h-',i1,' =',1p,1e11.3,
     &        ' Center x-',i1,' =',1p,1e11.3))

      end
