c$Id:$
      subroutine jdisp2(x,rlam,ebig,xa1,xb1,xah,xbh,
     &                  ra1,rb1,rah,rbh,j1,j2,pr,sr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]

c      Purpose: Form linear elastic spring control arrays
c         x(3,3)    - Reference coordinates of points & properties
c         rlam(3,*) - Rotation array for rigid bodies
c         ebig(3,3) - Orientation of slider axis
c         xa(3)     - Location of body 1 connection point
c         xb(3)     - Location of body 2 connection point
c         ra(3)     - Location of body 1 center of mass
c         rb(3)     - Location of body 2 center of mass
c         j1        - Number of rigid body 1
c         j2        - Number of rigid body 2

c      Outputs:
c         pr(*)     - Residual array for joint
c         sr(*,*)   - Tangent  array for joint
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ptdat6.h'

      integer   i,j,j1,j2
      real*8    l0, ln, l1, dlh, dln, dl1, fh, temp1
      real*8    ean(3), ea1(3), eah(3), ebn(3), eb1(3), ebh(3)
      real*8    x(3,3), rlam(9,6,*), ebig(3,3)
      real*8    xa1(3),xb1(3),xah(3),xbh(3), xan(3), xbn(3)
      real*8    rah(3),ra1(3),rbh(3),rb1(3), pr(3,7), sr(21,21)
      real*8    rcea1(3), rceb1(3), rceah(3), rcebh(3)

      save

c     Form updated basis vectors at t_n and at t_n+1 (e_a^3,e_b^3)
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


      do i = 1,3

c       Interpolate the base vectors to t_n+1/2
        eah(i) = 0.5d0*(ean(i) + ea1(i))
        ebh(i) = 0.5d0*(ebn(i) + eb1(i))

c       Interpolate the position vectors to t_n
        xan(i) = 2.d0*xah(i) - xa1(i)
        xbn(i) = 2.d0*xbh(i) - xb1(i)

      end do ! i

      l0 = 0.d0
      l1 = 0.d0
      ln = 0.d0

      do i=1,3
        l0 = l0 + (x(i,1) - x(i,2))*ebig(i,3)
        l1 = l1 + xa1(i)*ea1(i) - xb1(i)*eb1(i)
        ln = ln + xan(i)*ean(i) - xbn(i)*ebn(i)
      end do ! i

      dl1 = (l1 - l0)
      dln = (ln - l0)
      dlh = 0.5d0*(dl1 + dln)
      fh  = x(1,3)*dlh

      rceah(1) = rah(2)*eah(3) - rah(3)*eah(2)
      rceah(2) = rah(3)*eah(1) - rah(1)*eah(3)
      rceah(3) = rah(1)*eah(2) - rah(2)*eah(1)

      rcebh(1) = rbh(2)*ebh(3) - rbh(3)*ebh(2)
      rcebh(2) = rbh(3)*ebh(1) - rbh(1)*ebh(3)
      rcebh(3) = rbh(1)*ebh(2) - rbh(2)*ebh(1)

      do i = 1,3
        pr(i,1) = -fh*eah(i)
        pr(i,2) =  fh*rceah(i)
        pr(i,3) =  fh*ebh(i)
        pr(i,4) = -fh*rcebh(i)
      end do ! i

      rcea1(1) = ra1(2)*ea1(3) - ra1(3)*ea1(2)
      rcea1(2) = ra1(3)*ea1(1) - ra1(1)*ea1(3)
      rcea1(3) = ra1(1)*ea1(2) - ra1(2)*ea1(1)

      rceb1(1) = rb1(2)*eb1(3) - rb1(3)*eb1(2)
      rceb1(2) = rb1(3)*eb1(1) - rb1(1)*eb1(3)
      rceb1(3) = rb1(1)*eb1(2) - rb1(2)*eb1(1)

      do i = 1,3
        do j = 1,3
          sr(i,j  ) =  eah(i)*ea1(j)
          sr(i,j+3) = -eah(i)*rcea1(j)
          sr(i,j+6) = -eah(i)*eb1(j)
          sr(i,j+9) =  eah(i)*rceb1(j)

          sr(i+6,j  ) =  ebh(i)*ea1(j)
          sr(i+6,j+3) = -ebh(i)*rcea1(j)
          sr(i+6,j+6) = -ebh(i)*eb1(j)
          sr(i+6,j+9) =  ebh(i)*rceb1(j)
        end do ! j
      end do ! i

      sr(1,5) = sr(1,5) + dlh*ea1(3)
      sr(1,6) = sr(1,6) - dlh*ea1(2)
      sr(2,4) = sr(2,4) - dlh*ea1(3)
      sr(2,6) = sr(2,6) + dlh*ea1(1)
      sr(3,4) = sr(3,4) + dlh*ea1(2)
      sr(3,5) = sr(3,5) - dlh*ea1(1)

      sr(7,11) = sr(7,11) + dlh*eb1(3)
      sr(7,12) = sr(7,12) - dlh*eb1(2)
      sr(8,10) = sr(8,10) - dlh*eb1(3)
      sr(8,12) = sr(8,12) + dlh*eb1(1)
      sr(9,10) = sr(9,10) + dlh*eb1(2)
      sr(9,11) = sr(9,11) - dlh*eb1(1)

      temp1 = 0.d0
      do i = 1,3
        temp1 = temp1 + ea1(i)*rah(i)
      end do ! i

      do i = 1,3
        do j = 1,3
          sr(i+3,j  ) =  rceah(i)*ea1(j)
          sr(i+3,j+3) = -rceah(i)*rcea1(j) - dlh*ea1(i)*rah(j)
          sr(i+3,j+6) = -rceah(i)*eb1(j)
          sr(i+3,j+9) =  rceah(i)*rceb1(j)
        end do ! j
        sr(i+3,i+3) = sr(i+3,i+3) + dlh*temp1
      end do ! i

      sr(4,2) = sr(4,2) + dlh*eah(3)
      sr(4,3) = sr(4,3) - dlh*eah(2)
      sr(5,1) = sr(5,1) - dlh*eah(3)
      sr(5,3) = sr(5,3) + dlh*eah(1)
      sr(6,1) = sr(6,1) + dlh*eah(2)
      sr(6,2) = sr(6,2) - dlh*eah(1)

      temp1 = 0.d0
      do i = 1,3
        temp1 = temp1 + eb1(i)*rbh(i)
      end do ! i

      do i = 1,3
        do j = 1,3
          sr(i+9,j  ) =  rcebh(i)*ea1(j)
          sr(i+9,j+3) = -rcebh(i)*rcea1(j)
          sr(i+9,j+6) = -rcebh(i)*eb1(j)
          sr(i+9,j+9) =  rcebh(i)*rceb1(j) - dlh*eb1(i)*rbh(j)
        end do ! j
        sr(i+9,i+9) = sr(i+9,i+9) + dlh*temp1
      end do ! i

      sr(10,8) = sr(10,8) + dlh*ebh(3)
      sr(10,9) = sr(10,9) - dlh*ebh(2)
      sr(11,7) = sr(11,7) - dlh*ebh(3)
      sr(11,9) = sr(11,9) + dlh*ebh(1)
      sr(12,7) = sr(12,7) + dlh*ebh(2)
      sr(12,8) = sr(12,8) - dlh*ebh(1)

      do i = 1,3
        do j = 1,12
          sr(i  ,j) =  0.5d0*x(1,3)*sr(i  ,j)
          sr(i+3,j) = -0.5d0*x(1,3)*sr(i+3,j)
          sr(i+6,j) = -0.5d0*x(1,3)*sr(i+6,j)
          sr(i+9,j) =  0.5d0*x(1,3)*sr(i+9,j)
        end do ! j
      end do ! i

      epl(8) = epl(8) + 0.5d0*x(1,3)*(l1-l0)*(l1-l0)

      end
