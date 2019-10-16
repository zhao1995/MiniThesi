c$Id:$
      subroutine jdisp1(def,rlam,ebig,x1,x2,j1,j2,pr,sr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]

c      Purpose: Form FORCE control arrays

c      Inputs:
c         def(1)    - user defined control value
c         def(2)    - proportional load value applied to def(1)
c         rlam(3,*) - Array of rotation quantities for rigid bodies
c         ebig(3,3) - Orientation for revolute axis in reference frame
c         j1        - Number of rigid body 1
c         j2        - Number of rigid body 2

c      Outputs:
c         pr(*)     - Residual array for joint
c         sr(*,*)   - Tangent  array for joint
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'prld1.h'

      integer   i,j1,j2
      real*8    dforc, dist, ebig(3,3), rlam(9,6,*)
      real*8    ean(3), ea1(3), ebn(3), eb1(3)
      real*8    x1(3), x2(3), def(2), pr(3,7), sr(21,21)

      save

      dforc = def(1)*prldv(int(def(2)))
      dist  = sqrt((x1(1) - x2(1))*(x1(1) - x2(1))
     &           + (x1(2) - x2(2))*(x1(2) - x2(2))
     &           + (x1(3) - x2(3))*(x1(3) - x2(3)))

c     ea = e_3^a
c     eb = e_3^b

c     Form updated basis vectors for body 'j1'

      if(j1.gt.0) then
        call quavec(rlam(1,1,j1),ebig(1,3),ean)
        call quavec(rlam(1,3,j1),ebig(1,3),ea1)
      else
        do i = 1,3
          ean(i) = ebig(i,3)
          ea1(i) = ebig(i,3)
        end do ! i
      endif

c     Form updated basis vectors for body 'j2'

      if(j2.gt.0) then
        call quavec(rlam(1,1,j2),ebig(1,3),ebn)
        call quavec(rlam(1,3,j2),ebig(1,3),eb1)
      else
        do i = 1,3
          ebn(i) = ebig(i,3)
          eb1(i) = ebig(i,3)
        end do ! i
      endif

      sr(1,5) =  dforc*ea1(3)
      sr(1,6) = -dforc*ea1(2)
      sr(2,6) =  dforc*ea1(1)
      sr(2,4) = -dforc*ea1(3)
      sr(3,4) =  dforc*ea1(2)
      sr(3,5) = -dforc*ea1(1)

      sr(7,11) = -dforc*eb1(3)
      sr(7,12) =  dforc*eb1(2)
      sr(8,12) = -dforc*eb1(1)
      sr(8,10) =  dforc*eb1(3)
      sr(9,10) = -dforc*eb1(2)
      sr(9,11) =  dforc*eb1(1)

      do i = 1,3
        pr(i,1) = -dforc*ea1(i)
        pr(i,3) = +dforc*eb1(i)
      end do ! i

      end
