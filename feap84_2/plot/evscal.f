c$Id:$
      subroutine evscal(x,ct, es)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Scale factor for eigenvectors

c      Inputs:
c         x(ndm,*) - Nodal coordinate array
c         ct       - Scaling values

c      Output:
c         es     -  Scale of vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'

      integer    i,n
      real*8     x(ndm,numnp),ct,es, xmax,xmin

      save

      xmax = x(1,1)
      xmin = x(1,1)
      do n = 1,numnp
        do i = 1,ndm
          xmax = max(xmax,x(i,n))
          xmin = min(xmin,x(i,n))
        end do ! i
      end do ! n

      if(ct.eq.0.0d0) then
        es = 0.10d0*(xmax-xmin)
      else
        es = ct*(xmax-xmin)
      endif

      end
