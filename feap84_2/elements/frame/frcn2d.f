c$Id:$
      subroutine frcn2d(p,dt,st)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 2-d contour

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'strnum.h'

      integer   i,j
      real*8    dt(*),st(nen,*),p(3,*)

      save

      do i = 1,2

        dt(i) = dt(i) + 1.d0

c       Stress projections

        do j = 1,3
          st(i,j) = st(i,j) + p(j,i)
        end do ! j
      end do ! i

      iste = 3

      end
