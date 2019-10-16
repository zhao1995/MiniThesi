c$Id:$
      subroutine pcorner3d()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute deformation gradient and its inverse at tn+1

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'eldata.h'  ! nel
      include   'corner.h' ! ncorner,icorner(24,4)

      integer    jcorner4(2,12), jcorner8(2,24),  nn

      data       jcorner4/ 1,3, 1,2, 1,4,   2,1, 2,3, 2,4,
     &                     3,2, 3,1, 3,4,   4,1, 4,2, 4,3/

      data       jcorner8/ 1,3, 1,6, 1,8,   2,4, 2,5, 2,7,
     &                     3,1, 3,6, 3,8,   4,2, 4,5, 4,7,
     &                     5,2, 5,4, 5,7,   6,1, 6,3, 6,8,
     &                     7,2, 7,4, 7,5,   8,1, 8,3, 8,6/

c     data       jcorner8/ 1,4, 1,2, 1,5,   2,1, 2,3, 2,6,
c    &                     3,2, 3,4, 3,7,   4,3, 4,1, 4,8,
c    &                     5,8, 5,6, 5,1,   6,5, 6,7, 6,2,
c    &                     7,6, 7,8, 7,3,   8,7, 8,5, 8,4/

      if(nel.eq.4) then
        ncorner = 12
        ncorner = 0
        do nn = 1,ncorner
          icorner(nn, 1) = jcorner4(1,nn)
          icorner(nn, 2) = jcorner4(2,nn)
        end do ! nn

      else
        ncorner = 24
        do nn = 1,ncorner
          icorner(nn, 1) = jcorner8(1,nn)
          icorner(nn, 2) = jcorner8(2,nn)
        end do ! nn
      endif

      end
