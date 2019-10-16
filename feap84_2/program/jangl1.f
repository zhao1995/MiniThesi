c$Id:$
      subroutine jangl1(def,rlam,ebig,j1,j2,pr,sr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]

c      Purpose: Form TORQUE control arrays

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
      real*8    ebig(3,3), rlam(9,6,*), pr(3,7), sr(21,21)
      real*8    ean(3), ea1(3), eah(3), ebn(3), eb1(3), ebh(3)
      real*8    def(3), dtorqh, dtorq

      save

      dtorq = def(1)*prldv(int(def(2)))

      write(*,*)' JANGL1:dtorq = ',dtorq

c     ea = e_3^a
c     eb = e_3^b

c     Form updated basis vectors for body 'j1'

      if(j1.gt.0) then
        call quavec(rlam(1,1,j1),ebig(1,3),ean(1))
        call quavec(rlam(1,3,j1),ebig(1,3),ea1(1))
      else
        do i = 1,3
          ean(i) = ebig(i,3)
          ea1(i) = ebig(i,3)
        end do ! i
      endif

c     Form updated basis vectors for body 'j2'

      if(j2.gt.0) then
        call quavec(rlam(1,1,j2),ebig(1,3),ebn(1))
        call quavec(rlam(1,3,j2),ebig(1,3),eb1(1))
      else
        do i = 1,3
          ebn(i) = ebig(i,3)
          eb1(i) = ebig(i,3)
        end do ! i
      endif

c     Interpolate the base vectors to t_n+1/2

      do i = 1,3
        ebh(i) = 0.5d0*(ebn(i) + eb1(i))
        eah(i) = 0.5d0*(ean(i) + ea1(i))
      end do ! i

      dtorqh = 1.d0*dtorq

      sr(1,2) =  dtorqh*ea1(3)
      sr(1,3) = -dtorqh*ea1(2)
      sr(2,3) =  dtorqh*ea1(1)
      sr(2,2) = -dtorqh*ea1(3)
      sr(3,2) =  dtorqh*ea1(2)
      sr(3,2) = -dtorqh*ea1(1)

      sr(4,5) = -dtorqh*eb1(3)
      sr(4,6) =  dtorqh*eb1(2)
      sr(5,6) = -dtorqh*eb1(1)
      sr(5,4) =  dtorqh*eb1(3)
      sr(6,4) = -dtorqh*eb1(2)
      sr(6,5) =  dtorqh*eb1(1)

      do i = 1,3
        pr(i,1) = -dtorq*ea1(i)
        pr(i,2) =  dtorq*eb1(i)
      end do ! i

      end
