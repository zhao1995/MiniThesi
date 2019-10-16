c$Id:$
      subroutine pcorner2d()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute corner nodes

c     Inputs:
c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'eldata.h'  ! nel
      include   'corner.h' ! ncorner,icorner(24,4)
      include   'qudshp.h'

      integer    jcorner3(2,6),jcorner4(2,8),  nn

      data       jcorner3 / 1,3, 1,2,  2,1, 2,3, 3,2, 3,1 /
      data       jcorner4 / 1,4, 1,2,  2,1, 2,3, 3,2, 3,4, 4,3, 4,1 /

      if(nurbfl) return
      if(nel.eq.3 .or. nel.eq.6) then
        ncorner = 6
        do nn = 1,ncorner
          icorner(nn, 1) = jcorner3(1,nn)
          icorner(nn, 2) = jcorner3(2,nn)
        end do ! nn
      else
        ncorner = 8
        do nn = 1,ncorner
          icorner(nn, 1) = jcorner4(1,nn)
          icorner(nn, 2) = jcorner4(2,nn)
        end do ! nn
      endif

      end
