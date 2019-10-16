c$Id:$
      subroutine jangl2(stiff,rlam,ebig,j1,j2,pr,sr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]

c      Purpose: Form linear torsional spring control arrays.

c      Inputs:
c         stiff     - Value of lagrange multiplier force on joint
c         rlam(3,*) - Array of rotation quantities for rigid bodies
c         ebig(3,3) - Orientation for revolute axis in reference frame
c         j1        - Number of rigid body 1
c         j2        - Number of rigid body 2

c      Outputs:
c         pr(*)     - Residual array for joint
c         sr(*,*)   - Tangent  array for joint
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'aceang.h'
      include  'ptdat6.h'

      integer   i,j,j1,j2
      real*8    stiff, ebig(3,3), rlam(9,6,*), pr(3,7), sr(21,21)
      real*8    ean(3), ea1(3), eah(3), ebn(3), eb1(3), ebh(3)
      real*8    mcona(3,3), mconb(3,3), cona(3), conb(3)
      real*8    mh, danga, dangb

      save

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

c     Interpolate the base vectors to t_n+1/2

      do i = 1,3
        eah(i) = 0.5d0*(ean(i) + ea1(i))
        ebh(i) = 0.5d0*(ebn(i) + eb1(i))
      end do ! i

      danga = 0.d0
      dangb = 0.d0
      do i = 1,3
        danga = danga + rlam(i,5,j1)*ea1(i)
        dangb = dangb + rlam(i,5,j2)*eb1(i)
      end do ! i

      ang(2) = ang(1) + dangb - danga
      mh = ang(2) + ang(1)

      do i = 1,3
        mcona(i,i) = 1.d0
        mconb(i,i) = 1.d0
        cona(i) = 0.d0
        conb(i) = 0.d0
      end do ! i

      mcona(1,2) = rlam(3,5,j1)
      mcona(1,3) = rlam(2,5,j1)
      mcona(2,1) = rlam(3,5,j1)
      mcona(2,3) = rlam(1,5,j1)
      mcona(3,1) = rlam(2,5,j1)
      mcona(3,2) = rlam(1,5,j1)

      mconb(1,2) = rlam(3,5,j2)
      mconb(1,3) = rlam(2,5,j2)
      mconb(2,1) = rlam(3,5,j2)
      mconb(2,3) = rlam(1,5,j2)
      mconb(3,1) = rlam(2,5,j2)
      mconb(3,2) = rlam(1,5,j2)

      do i = 1,3
        do j = 1,3
         cona(i) = cona(i) + mcona(i,j)*ea1(j)
         conb(i) = conb(i) + mconb(i,j)*eb1(j)
        end do ! j
      end do ! i

      do i = 1,3
        do j = 1,3
          sr(i  ,j  ) = -eah(i)*cona(j)
          sr(i  ,j+3) =  eah(i)*conb(j)
          sr(i+3,j  ) = -ebh(i)*cona(j)
          sr(i+3,j+3) =  ebh(i)*conb(j)
        end do ! j
      end do ! i

      sr(1,2) = sr(1,2) + mh*ea1(3)
      sr(1,3) = sr(1,3) - mh*ea1(2)
      sr(2,1) = sr(2,1) - mh*ea1(3)
      sr(2,3) = sr(2,3) + mh*ea1(1)
      sr(3,1) = sr(3,1) + mh*ea1(2)
      sr(3,2) = sr(3,2) - mh*ea1(1)

      sr(4,5) = sr(4,5) + mh*eb1(3)
      sr(4,6) = sr(4,6) - mh*eb1(2)
      sr(5,4) = sr(5,4) - mh*eb1(3)
      sr(5,6) = sr(5,6) + mh*eb1(1)
      sr(6,4) = sr(6,4) + mh*eb1(2)
      sr(6,5) = sr(6,5) - mh*eb1(1)

      do i = 1,3
        do j = 1,6
          sr(i  ,j) = -0.5d0*stiff*sr(i  ,j)
          sr(i+3,j) =  0.5d0*stiff*sr(i+3,j)
        end do ! j
      end do ! i

      do i = 1,3
        pr(i,1) =  0.5d0*stiff*mh*eah(i)
        pr(i,2) = -0.5d0*stiff*mh*ebh(i)
      end do ! i

      epl(7) = epl(7) + 0.5d0*stiff*ang(2)*ang(2)

      end
