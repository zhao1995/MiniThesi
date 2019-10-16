c$Id:$
      subroutine jangl0(lam,angle,def,rlam,ebig,j1,j2,pr,sr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved
c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]

c      Purpose: Form ANGULAR control arrays

c      Inputs:
c         lam       - Value of lagrange multiplier force on joint
c         angle     - Angular offset between body frames (in rads)
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

      include  'debugs.h'
      include  'pconstant.h'
      include  'prld1.h'

      integer   i,j,j1,j2
      real*8    lam,lamh, angle,dangle, ebig(3,3), rlam(9,6,*)
      real*8    ean(3,2),ea1(3,2),eah(3,2), ebn(3),eb1(3),ebh(3)
      real*8    pr(3,7), sr(21,21), d(8), c(8,3), def(3)
      real*8    zetaa(3,3), zetab(3,3), gam1(3), gamh(3)

      save

      dangle = def(1)*prldv(nint(def(2)))

      do while(abs(dangle - angle).gt.pi)
        if(dangle.gt.0.d0) then
          angle = angle + 2.d0*pi
        elseif(dangle.lt.0.d0) then
          angle = angle - 2.d0*pi
        end if
      end do ! while

      if(debug) then
        write(*,*)'JTHETA:lam = ',lam,' and dangle = ',dangle
      endif

c     Form updated basis vectors for body 'j1'

      if(j1.gt.0) then
        call quavec(rlam(1,1,j1),ebig(1,1),ean(1,1))
        call quavec(rlam(1,1,j1),ebig(1,2),ean(1,2))
        call quavec(rlam(1,3,j1),ebig(1,1),ea1(1,1))
        call quavec(rlam(1,3,j1),ebig(1,2),ea1(1,2))
      else
        do i = 1,3
          ean(i,1) = ebig(i,1)
          ean(i,2) = ebig(i,2)
          ea1(i,1) = ebig(i,1)
          ea1(i,2) = ebig(i,2)
        end do ! i
      endif

c     Form updated basis vectors for body 'j2'

      if(j2.gt.0) then
        call quavec(rlam(1,1,j2),ebig(1,1),ebn(1))
        call quavec(rlam(1,3,j2),ebig(1,1),eb1(1))
      else
        do i = 1,3
          ebn(i) = ebig(i,1)
          eb1(i) = ebig(i,1)
        end do ! i
      endif

c     Interpolate the base vectors to t_n+1/2

      do i = 1,3
        ebh(i) = 0.5d0*(ebn(i) + eb1(i))
        do j = 1,2
          eah(i,j) = 0.5d0*(ean(i,j) + ea1(i,j))
        end do ! j
      end do ! i

      call pzero(d,8)
      do i = 1,3
        d(1)  = d(1) + ebh(i)*eah(i,1)
        d(2)  = d(2) + ebh(i)*eah(i,2)
        d(3)  = d(3) + eb1(i)*ea1(i,1)
        d(4)  = d(4) + eb1(i)*ea1(i,2)
        d(5)  = d(5) + ebh(i)*ea1(i,1)
        d(6)  = d(6) + ebh(i)*ea1(i,2)
        d(7)  = d(7) + eb1(i)*eah(i,1)
        d(8)  = d(8) + eb1(i)*eah(i,2)
      end do ! i

      c(1,1) = ebh(2)*eah(3,1) - ebh(3)*eah(2,1)
      c(1,2) = ebh(3)*eah(1,1) - ebh(1)*eah(3,1)
      c(1,3) = ebh(1)*eah(2,1) - ebh(2)*eah(1,1)

      c(2,1) = ebh(2)*eah(3,2) - ebh(3)*eah(2,2)
      c(2,2) = ebh(3)*eah(1,2) - ebh(1)*eah(3,2)
      c(2,3) = ebh(1)*eah(2,2) - ebh(2)*eah(1,2)

      c(3,1) = eb1(2)*ea1(3,1) - eb1(3)*ea1(2,1)
      c(3,2) = eb1(3)*ea1(1,1) - eb1(1)*ea1(3,1)
      c(3,3) = eb1(1)*ea1(2,1) - eb1(2)*ea1(1,1)

      c(4,1) = eb1(2)*ea1(3,2) - eb1(3)*ea1(2,2)
      c(4,2) = eb1(3)*ea1(1,2) - eb1(1)*ea1(3,2)
      c(4,3) = eb1(1)*ea1(2,2) - eb1(2)*ea1(1,2)

      c(5,1) = ebh(2)*ea1(3,1) - ebh(3)*ea1(2,1)
      c(5,2) = ebh(3)*ea1(1,1) - ebh(1)*ea1(3,1)
      c(5,3) = ebh(1)*ea1(2,1) - ebh(2)*ea1(1,1)

      c(6,1) = ebh(2)*ea1(3,2) - ebh(3)*ea1(2,2)
      c(6,2) = ebh(3)*ea1(1,2) - ebh(1)*ea1(3,2)
      c(6,3) = ebh(1)*ea1(2,2) - ebh(2)*ea1(1,2)

      c(7,1) = eb1(2)*eah(3,1) - eb1(3)*eah(2,1)
      c(7,2) = eb1(3)*eah(1,1) - eb1(1)*eah(3,1)
      c(7,3) = eb1(1)*eah(2,1) - eb1(2)*eah(1,1)

      c(8,1) = eb1(2)*eah(3,2) - eb1(3)*eah(2,2)
      c(8,2) = eb1(3)*eah(1,2) - eb1(1)*eah(3,2)
      c(8,3) = eb1(1)*eah(2,2) - eb1(2)*eah(1,2)

      do i = 1,3
        gam1(i) = d(4)*c(3,i) - d(3)*c(4,i)
        gamh(i) = d(2)*c(1,i) - d(1)*c(2,i)
      end do ! i

      do i = 1,3
        do j = 1,3
          zetaa(i,j) = - c(1,i)*c(6,j) + c(2,i)*c(5,j)
     &                 - d(2)*ea1(i,1)*ebh(j)
     &                 + d(1)*ea1(i,2)*ebh(j)
          zetab(i,j) =   c(1,i)*c(8,j) - c(2,i)*c(7,j)
     &                 + d(2)*eb1(i)*eah(j,1)
     &                 - d(1)*eb1(i)*eah(j,2)
        end do ! j
        zetaa(i,i) = zetaa(i,i) + d(2)*d(5) - d(1)*d(6)
        zetab(i,i) = zetab(i,i) - d(2)*d(7) + d(1)*d(8)
      end do ! i

      lamh = 0.5d0*lam

      do i = 1,3
        do j = 1,3
          sr(i  ,j  ) = sr(i  ,j  ) + lamh*zetaa(i,j)
          sr(i  ,j+3) = sr(i  ,j+3) - lamh*zetaa(i,j)
          sr(i+3,j  ) = sr(i+3,j  ) + lamh*zetab(i,j)
          sr(i+3,j+3) = sr(i+3,j+3) - lamh*zetab(i,j)
        end do ! j
        sr(i  ,7  ) = sr(i  ,7  ) + gamh(i)
        sr(i+3,7  ) = sr(i+3,7  ) - gamh(i)
        sr(7  ,i  ) = sr(7  ,i  ) + gam1(i)
        sr(7  ,i+3) = sr(7  ,i+3) - gam1(i)
      end do ! i

      do i = 1,3
        pr(i,1) = pr(i,1) - lam*gamh(i)
        pr(i,2) = pr(i,2) + lam*gamh(i)
      end do ! i
      pr(1,3) = pr(1,3) - (angle - dangle)

      end
