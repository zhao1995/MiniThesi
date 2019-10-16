c$Id:$
      subroutine pltelnd(x,ix,ndm,nel,num)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot location of nodal points in mesh

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         ndm       - Dimension of x array
c         nel       - Number of nodes on element
c         num       - Plot numbers if true

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata1.h'

      logical   num
      integer   ix(*),    ndm, nel, n
      real*8    x(ndm,*), x3, dx1

      save

c     Plot locations of all nodes

      dx1 = 0.002d0/scale
      x3  = 0.0d0
      do n = 1,nel
        if(ndm.ge.3) x3 = x(3,n)
        call plotl(x(1,n)-dx1 , x(2,n)+dx1 , x3, 3)
        call plotl(x(1,n)-dx1 , x(2,n)-dx1 , x3, 2)
        call plotl(x(1,n)+dx1 , x(2,n)-dx1 , x3, 2)
        call plotl(x(1,n)+dx1 , x(2,n)+dx1 , x3, 2)
        call plotl(x(1,n)-dx1 , x(2,n)+dx1 , x3, 2)
        call plotl(x(1,n)+2.d0*dx1, x(2,n)+2.d0*dx1, x3+1.d0*dx1, 3)
        if(num) call plabl(ix(n))
      end do ! n

      end
